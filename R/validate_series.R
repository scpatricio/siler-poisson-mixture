# Numerical checks for the Siler remaining-life-expectancy formulas.
#
# The script compares adaptive quadrature of conditional Siler survival with
# the Poisson-weighted generalized-exponential-integral representation.

scaled_expint = function(nu, z, rel.tol = 1e-12) {
  stopifnot(nu > 0, z > 0)
  if (!requireNamespace("expint", quietly = TRUE)) {
    stop("Install package 'expint' before running this validation.")
  }
  value = exp(z) * z^(nu - 1) * expint::gammainc(1 - nu, z)
  if (is.finite(value) && value > 0) {
    return(value)
  }
  integrate(
    function(s) exp(-z * expm1(s) - (nu - 1) * s),
    lower = 0,
    upper = Inf,
    rel.tol = rel.tol,
    subdivisions = 2000L,
    stop.on.error = TRUE
  )$value
}

siler_rle_quadrature = function(x, a1, b1, a2, a3, b3,
                                 rel.tol = 1e-11) {
  stopifnot(x >= 0, a1 >= 0, b1 > 0, a2 >= 0, a3 > 0, b3 > 0)
  A = a1 * exp(-b1 * x) / b1
  G = a3 * exp(b3 * x) / b3

  conditional_survival = function(t) {
    bt = b3 * t
    old_age_increment = rep(Inf, length(bt))
    finite_increment = bt <= 700
    old_age_increment[finite_increment] = expm1(bt[finite_increment])
    log_p = -A * (-expm1(-b1 * t)) - a2 * t - G * old_age_increment
    exp(log_p)
  }

  integrate(
    conditional_survival,
    lower = 0,
    upper = Inf,
    rel.tol = rel.tol,
    subdivisions = 2000L,
    stop.on.error = TRUE
  )$value
}

siler_rle_series = function(x, a1, b1, a2, a3, b3,
                             absolute.tol = 1e-12) {
  stopifnot(x >= 0, a1 >= 0, b1 > 0, a2 >= 0, a3 > 0, b3 > 0)
  A = a1 * exp(-b1 * x) / b1
  G = a3 * exp(b3 * x) / b3
  r = a2 / b3
  d = b1 / b3

  # The exact representation can be written as an expectation over
  # N ~ Poisson(A).  This form avoids the potentially unstable factor exp(-A).
  n_max = 0L
  repeat {
    tail_probability = ppois(n_max, A, lower.tail = FALSE)
    absolute_error_bound = tail_probability / (a2 + (n_max + 1) * b1)
    if (absolute_error_bound <= absolute.tol) {
      break
    }
    n_max = n_max + 1L
  }
  n = 0:n_max
  weights = dpois(n, lambda = A)
  values = vapply(
    n,
    function(j) scaled_expint(1 + r + j * d, G),
    numeric(1)
  )

  estimate = sum(weights * values) / b3
  gm_upper = scaled_expint(1 + r, G) / b3
  absolute_error_bound = ppois(n_max, A, lower.tail = FALSE) /
    (a2 + (n_max + 1) * b1)

  list(
    estimate = estimate,
    terms = n_max + 1L,
    absolute_error_bound = absolute_error_bound,
    gm_lower = exp(-A) * gm_upper,
    gm_upper = gm_upper,
    A = A,
    G = G
  )
}

parameter_sets = list(
  low_background = c(a1 = 0.15, b1 = 1.20, a2 = 0.002,
                     a3 = 0.00002, b3 = 0.10),
  elevated_early = c(a1 = 0.40, b1 = 0.30, a2 = 0.008,
                     a3 = 0.00008, b3 = 0.085),
  intermediate = c(a1 = 0.25, b1 = 0.50, a2 = 0.010,
                   a3 = 0.00010, b3 = 0.080)
)

ages = c(0, 30, 65)

rows = do.call(rbind, lapply(names(parameter_sets), function(label) {
  par = parameter_sets[[label]]
  do.call(rbind, lapply(ages, function(x) {
    series = do.call(siler_rle_series, c(list(x = x), as.list(par)))
    quadrature = do.call(siler_rle_quadrature, c(list(x = x), as.list(par)))
    data.frame(
      scenario = label,
      age = x,
      A_x = series$A,
      quadrature = quadrature,
      series = series$estimate,
      absolute_difference = abs(series$estimate - quadrature),
      terms = series$terms,
      bound = series$absolute_error_bound,
      gm_lower = series$gm_lower,
      gm_upper = series$gm_upper,
      row.names = NULL
    )
  }))
}))

print(rows, digits = 10)

stopifnot(all(rows$series >= rows$gm_lower - 1e-9))
stopifnot(all(rows$series <= rows$gm_upper + 1e-9))
stopifnot(all(rows$absolute_difference < 1e-8))

