# Core functions for Siler mortality, Poisson maximum likelihood estimation,
# remaining life expectancy, and diagnostics.

required_columns = c("year", "age", "deaths", "exposure")

check_analysis_data = function(data) {
  missing_columns = setdiff(required_columns, names(data))
  if (length(missing_columns) > 0L) {
    stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
  }
  if (anyDuplicated(data[c("year", "age")])) {
    stop("The analysis data contain duplicated year-age keys.")
  }
  if (any(!is.finite(data$deaths)) || any(data$deaths < 0)) {
    stop("Death counts must be finite and nonnegative.")
  }
  if (any(!is.finite(data$exposure)) || any(data$exposure <= 0)) {
    stop("Exposures must be finite and strictly positive.")
  }
  invisible(TRUE)
}

theta_names = c("a1", "b1", "a2", "a3", "b3")

theta_from_eta = function(eta) {
  theta = exp(eta)
  names(theta) = theta_names
  theta
}

eta_from_theta = function(theta) {
  theta = theta[theta_names]
  if (any(!is.finite(theta)) || any(theta <= 0)) {
    stop("All Siler parameters must be finite and strictly positive.")
  }
  eta = log(theta)
  names(eta) = theta_names
  eta
}

log_sum_exp_three = function(x1, x2, x3) {
  m = pmax(x1, x2, x3)
  m + log(exp(x1 - m) + exp(x2 - m) + exp(x3 - m))
}

siler_log_hazard = function(age, eta) {
  theta = theta_from_eta(eta)
  log_juvenile = eta["a1"] - theta["b1"] * age
  log_background = rep(eta["a2"], length(age))
  log_senescent = eta["a3"] + theta["b3"] * age
  log_sum_exp_three(log_juvenile, log_background, log_senescent)
}

siler_hazard_components = function(age, theta) {
  data.frame(
    age = age,
    juvenile = theta["a1"] * exp(-theta["b1"] * age),
    background = rep(theta["a2"], length(age)),
    senescent = theta["a3"] * exp(theta["b3"] * age),
    check.names = FALSE
  )
}

siler_hazard = function(age, theta) {
  components = siler_hazard_components(age, theta)
  components$juvenile + components$background + components$senescent
}

siler_nll = function(eta, age, deaths, exposure) {
  transformed = exp(eta)
  if (any(!is.finite(transformed)) || any(transformed == 0)) {
    return(.Machine$double.xmax / 100)
  }
  log_mu = siler_log_hazard(age, eta)
  log_expected = log(exposure) + log_mu
  if (any(!is.finite(log_expected)) || any(log_expected > 700)) {
    return(.Machine$double.xmax / 100)
  }
  expected = exp(log_expected)
  value = sum(expected - deaths * log_mu)
  if (!is.finite(value)) {
    value = .Machine$double.xmax / 100
  }
  value
}

safe_log_rate = function(rate, floor_value = 1e-12) {
  log(pmax(rate, floor_value))
}

initial_siler_theta = function(data) {
  age = data$age
  rate = data$deaths / data$exposure
  positive_rate = rate[is.finite(rate) & rate > 0]
  rate_floor = if (length(positive_rate) > 0L) min(positive_rate) / 10 else 1e-10

  middle = age >= 15 & age <= 35 & is.finite(rate)
  a2 = median(rate[middle], na.rm = TRUE)
  if (!is.finite(a2) || a2 <= 0) {
    a2 = 0.001
  }

  old = age >= 60 & age <= 90 & rate > a2
  if (sum(old) >= 8L) {
    old_fit = lm(safe_log_rate(rate[old] - a2, rate_floor) ~ age[old])
    b3 = max(0.01, min(0.25, unname(coef(old_fit)[2])))
    a3 = exp(unname(coef(old_fit)[1]))
  } else {
    b3 = 0.09
    a3 = 1e-5
  }

  early_residual = rate - a2 - a3 * exp(b3 * age)
  early = age >= 0 & age <= 10 & early_residual > 0
  if (sum(early) >= 5L) {
    early_fit = lm(safe_log_rate(early_residual[early], rate_floor) ~ age[early])
    b1 = max(0.05, min(5.0, -unname(coef(early_fit)[2])))
    a1 = exp(unname(coef(early_fit)[1]))
  } else {
    b1 = 1.0
    a1 = max(rate[age == 0] - a2 - a3, 0.01)
  }

  theta = c(a1 = a1, b1 = b1, a2 = a2, a3 = a3, b3 = b3)
  theta = pmax(theta, c(a1 = 1e-8, b1 = 0.01, a2 = 1e-8, a3 = 1e-12, b3 = 0.001))
  theta
}

make_multistart_values = function(initial_eta, number_of_starts, seed) {
  set.seed(seed)
  starts = vector("list", number_of_starts)
  starts[[1L]] = initial_eta
  if (number_of_starts > 1L) {
    for (i in 2:number_of_starts) {
      perturbation = rnorm(length(initial_eta), mean = 0, sd = 0.75)
      starts[[i]] = initial_eta + perturbation
    }
  }
  starts
}

fit_siler_poisson = function(data, number_of_starts = 16L,
                              max_iterations = 5000L, seed = 1L) {
  check_analysis_data(data)
  initial_theta = initial_siler_theta(data)
  initial_eta = eta_from_theta(initial_theta)
  starts = make_multistart_values(
    initial_eta = initial_eta,
    number_of_starts = number_of_starts,
    seed = seed
  )

  fits = lapply(starts, function(start) {
    tryCatch(
      optim(
        par = start,
        fn = siler_nll,
        age = data$age,
        deaths = data$deaths,
        exposure = data$exposure,
        method = "Nelder-Mead",
        control = list(
          maxit = max_iterations,
          reltol = 1e-10
        )
      ),
      error = function(error) NULL
    )
  })

  valid = vapply(
    fits,
    function(fit) !is.null(fit) && is.finite(fit$value),
    logical(1)
  )
  if (!any(valid)) {
    stop("All Siler optimization starts failed.")
  }
  fits = fits[valid]
  converged = vapply(fits, function(fit) fit$convergence == 0L, logical(1))
  candidate_fits = if (any(converged)) fits[converged] else fits
  best = candidate_fits[[which.min(vapply(
    candidate_fits,
    function(fit) fit$value,
    numeric(1)
  ))]]
  eta_hat = best$par
  theta_hat = theta_from_eta(eta_hat)

  list(
    theta = theta_hat,
    eta = eta_hat,
    convergence = best$convergence,
    message = if (is.null(best$message)) "" else best$message,
    nll_kernel = best$value,
    starts_attempted = number_of_starts,
    starts_successful = length(fits),
    observations = nrow(data),
    total_deaths = sum(data$deaths),
    total_exposure = sum(data$exposure)
  )
}

scaled_generalized_expint_fallback = function(
    nu, z, relative_tolerance = 1e-11) {
  if (!is.finite(nu) || nu <= 0 || !is.finite(z) || z <= 0) {
    stop("The generalized exponential integral requires nu > 0 and z > 0.")
  }
  integrand = function(s) {
    value = numeric(length(s))
    finite_exponential = s < 700
    value[finite_exponential] = exp(
      -z * expm1(s[finite_exponential]) -
        (nu - 1) * s[finite_exponential]
    )
    value
  }
  integrate(
    integrand,
    lower = 0,
    upper = Inf,
    rel.tol = relative_tolerance,
    subdivisions = 2000L,
    stop.on.error = TRUE
  )$value
}

scaled_generalized_expint = function(nu, z, relative_tolerance = 1e-11) {
  if (!is.finite(nu) || nu <= 0 || !is.finite(z) || z <= 0) {
    stop("The generalized exponential integral requires nu > 0 and z > 0.")
  }
  if (!requireNamespace("expint", quietly = TRUE)) {
    stop(
      "Package 'expint' is required. Install it with install.packages('expint')."
    )
  }
  incomplete_gamma = suppressWarnings(expint::gammainc(1 - nu, z))
  unscaled_value = z^(nu - 1) * incomplete_gamma
  scaled_value = exp(z) * unscaled_value
  if (is.finite(scaled_value) && scaled_value > 0) {
    return(scaled_value)
  }
  scaled_generalized_expint_fallback(
    nu = nu,
    z = z,
    relative_tolerance = relative_tolerance
  )
}

siler_remaining_life = function(age, theta, absolute_tolerance = 1e-12,
                                 relative_tolerance = 1e-11) {
  A = theta["a1"] * exp(-theta["b1"] * age) / theta["b1"]
  G = theta["a3"] * exp(theta["b3"] * age) / theta["b3"]
  r = theta["a2"] / theta["b3"]
  d = theta["b1"] / theta["b3"]
  if (!is.finite(absolute_tolerance) || absolute_tolerance <= 0) {
    stop("absolute_tolerance must be positive and finite.")
  }
  n_max = 0L
  repeat {
    tail_probability = ppois(n_max, lambda = A, lower.tail = FALSE)
    truncation_bound = tail_probability /
      (theta["a2"] + (n_max + 1) * theta["b1"])
    if (truncation_bound <= absolute_tolerance) {
      break
    }
    n_max = n_max + 1L
  }
  n = 0:n_max
  weights = dpois(n, lambda = A)
  scaled_values = vapply(
    n,
    function(index) {
      scaled_generalized_expint(
        nu = 1 + r + index * d,
        z = G,
        relative_tolerance = relative_tolerance
      )
    },
    numeric(1)
  )
  gm_upper = scaled_generalized_expint(
    nu = 1 + r,
    z = G,
    relative_tolerance = relative_tolerance
  ) / theta["b3"]
  estimate = sum(weights * scaled_values) / theta["b3"]
  tail_probability = ppois(n_max, lambda = A, lower.tail = FALSE)
  truncation_bound = tail_probability /
    (theta["a2"] + (n_max + 1) * theta["b1"])

  result = c(
    estimate,
    exp(-A) * gm_upper,
    gm_upper,
    truncation_bound,
    n_max + 1L,
    A,
    G
  )
  names(result) = c(
    "exact",
    "gm_lower",
    "gm_upper",
    "truncation_bound",
    "terms",
    "A",
    "G"
  )
  result
}

siler_remaining_life_partial = function(age, theta, maximum_indices = 0:3,
                                         relative_tolerance = 1e-11) {
  maximum_indices = as.integer(maximum_indices)
  if (length(maximum_indices) == 0L || any(maximum_indices < 0L)) {
    stop("maximum_indices must contain nonnegative integers.")
  }
  A = theta["a1"] * exp(-theta["b1"] * age) / theta["b1"]
  G = theta["a3"] * exp(theta["b3"] * age) / theta["b3"]
  r = theta["a2"] / theta["b3"]
  d = theta["b1"] / theta["b3"]
  n = 0:max(maximum_indices)
  terms = dpois(n, lambda = A) * vapply(
    n,
    function(index) {
      scaled_generalized_expint(
        nu = 1 + r + index * d,
        z = G,
        relative_tolerance = relative_tolerance
      ) / theta["b3"]
    },
    numeric(1)
  )
  partial_sums = cumsum(terms)
  values = partial_sums[maximum_indices + 1L]
  names(values) = paste0("m", maximum_indices)
  values
}

siler_remaining_life_truncation_difference = function(
    age, theta, maximum_indices = 0:4, tail_tolerance = 1e-14,
    relative_tolerance = 1e-11) {
  maximum_indices = as.integer(maximum_indices)
  if (length(maximum_indices) == 0L || any(maximum_indices < 0L)) {
    stop("maximum_indices must contain nonnegative integers.")
  }
  A = theta["a1"] * exp(-theta["b1"] * age) / theta["b1"]
  G = theta["a3"] * exp(theta["b3"] * age) / theta["b3"]
  r = theta["a2"] / theta["b3"]
  d = theta["b1"] / theta["b3"]
  probability_limit = max(
    0L,
    as.integer(qpois(1 - tail_tolerance, lambda = A))
  )
  # Fixed-m errors require terms beyond the largest retained index even when
  # their total probability is already below the global series tolerance.
  # The padding preserves these small positive tails for logarithmic plots.
  tail_padding = 20L
  n_max = max(
    probability_limit,
    max(maximum_indices) + tail_padding
  )
  n = 0:n_max
  terms = dpois(n, lambda = A) * vapply(
    n,
    function(index) {
      scaled_generalized_expint(
        nu = 1 + r + index * d,
        z = G,
        relative_tolerance = relative_tolerance
      ) / theta["b3"]
    },
    numeric(1)
  )
  absolute_differences = vapply(
    maximum_indices,
    function(index) sum(terms[n > index]),
    numeric(1)
  )
  names(absolute_differences) = paste0("m", maximum_indices)
  absolute_differences
}

siler_remaining_life_quadrature = function(age, theta,
                                            relative_tolerance = 1e-10) {
  A = theta["a1"] * exp(-theta["b1"] * age) / theta["b1"]
  G = theta["a3"] * exp(theta["b3"] * age) / theta["b3"]
  conditional_survival = function(t) {
    scaled_time = theta["b3"] * t
    old_age_increment = rep(Inf, length(scaled_time))
    finite_increment = scaled_time <= 700
    old_age_increment[finite_increment] = expm1(scaled_time[finite_increment])
    log_value = -A * (-expm1(-theta["b1"] * t)) -
      theta["a2"] * t - G * old_age_increment
    exp(log_value)
  }
  integrate(
    conditional_survival,
    lower = 0,
    upper = Inf,
    rel.tol = relative_tolerance,
    subdivisions = 2000L,
    stop.on.error = TRUE
  )$value
}

format_fit_summary = function(fit, year, sex) {
  data.frame(
    year = year,
    sex = sex,
    a1 = fit$theta["a1"],
    b1 = fit$theta["b1"],
    a2 = fit$theta["a2"],
    a3 = fit$theta["a3"],
    b3 = fit$theta["b3"],
    convergence = fit$convergence,
    observations = fit$observations,
    total_deaths = fit$total_deaths,
    total_exposure = fit$total_exposure,
    starts_successful = fit$starts_successful,
    starts_attempted = fit$starts_attempted,
    row.names = NULL
  )
}

