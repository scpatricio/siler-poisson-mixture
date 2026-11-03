# Fit the Siler model to Swedish female mortality using a Poisson likelihood.
# This script reads Human Mortality Database death counts and exposures,
# fits the selected years in parallel, and saves analysis-ready results.

r_directory = file.path(getwd(), "R")
source(file.path(r_directory, "config.R"))
source(file.path(r_directory, "siler_functions.R"))

read_hmd_matrix = function(path, value_name, sex) {
  raw = read.table(
    file = path,
    skip = 2L,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = "."
  )
  required = c("Year", "Age", sex)
  if (!all(required %in% names(raw))) {
    stop("The file ", path, " does not contain: ", paste(required, collapse = ", "))
  }
  age_text = as.character(raw$Age)
  age = suppressWarnings(as.integer(age_text))
  data.frame(
    year = as.integer(raw$Year),
    age_text = age_text,
    age = age,
    value = as.numeric(raw[[sex]]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) |>
    setNames(c("year", "age_text", "age", value_name))
}

read_hmd_life_table = function(path) {
  raw = read.table(
    file = path,
    skip = 2L,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = "."
  )
  data.frame(
    year = as.integer(raw$Year),
    age = suppressWarnings(as.integer(as.character(raw$Age))),
    hmd_mx = as.numeric(raw$mx),
    hmd_ex = as.numeric(raw$ex),
    stringsAsFactors = FALSE
  )
}

sex = analysis_config$sex
if (!identical(sex, "Female")) {
  stop("This analysis is restricted to Swedish female mortality data.")
}

death_file = project_path("data", "hmd", "Deaths_1x1.txt")
exposure_file = project_path("data", "hmd", "Exposures_1x1.txt")
life_table_file = project_path("data", "hmd", "fltper_1x1.txt")

deaths = read_hmd_matrix(death_file, "deaths", sex)
exposures = read_hmd_matrix(exposure_file, "exposure", sex)

if (!identical(deaths[c("year", "age_text")], exposures[c("year", "age_text")])) {
  stop("Death and exposure files do not have identical year-age keys.")
}

analysis_data = merge(
  deaths[c("year", "age_text", "age", "deaths")],
  exposures[c("year", "age_text", "exposure")],
  by = c("year", "age_text"),
  all = FALSE,
  sort = TRUE
)
analysis_data = analysis_data[
  analysis_data$year %in% analysis_config$years &
    !is.na(analysis_data$age) &
    analysis_data$age >= analysis_config$minimum_age &
    analysis_data$age <= analysis_config$maximum_age &
    analysis_data$exposure > 0,
  c("year", "age", "deaths", "exposure")
]

if (!all(analysis_config$years %in% unique(analysis_data$year))) {
  missing_years = setdiff(analysis_config$years, unique(analysis_data$year))
  stop("Requested years are absent after filtering: ", paste(missing_years, collapse = ", "))
}
check_analysis_data(analysis_data)

life_table = read_hmd_life_table(life_table_file)
life_table = life_table[
  life_table$year %in% analysis_config$years &
    !is.na(life_table$age) &
    life_table$age %in% analysis_config$life_expectancy_ages,
  ,
  drop = FALSE
]

fit_one_year = function(year) {
  tryCatch(
    {
      year_data = analysis_data[analysis_data$year == year, , drop = FALSE]
      fit = fit_siler_poisson(
        data = year_data,
        number_of_starts = analysis_config$number_of_starts,
        max_iterations = analysis_config$optimizer_max_iterations,
        seed = analysis_config$random_seed + year
      )
      list(year = year, data = year_data, fit = fit, error = NULL)
    },
    error = function(error) {
      list(year = year, data = NULL, fit = NULL, error = conditionMessage(error))
    }
  )
}

years = analysis_config$years
number_of_cores = min(length(years), analysis_config$parallel_cores)
if (.Platform$OS.type == "unix" && number_of_cores > 1L) {
  fitted = parallel::mclapply(
    years,
    fit_one_year,
    mc.cores = number_of_cores,
    mc.preschedule = TRUE,
    mc.set.seed = TRUE
  )
} else {
  fitted = lapply(years, fit_one_year)
}
names(fitted) = as.character(years)

failed = vapply(fitted, function(item) !is.null(item$error), logical(1))
if (any(failed)) {
  failure_messages = vapply(
    fitted[failed],
    function(item) paste0(item$year, ": ", item$error),
    character(1)
  )
  stop(
    "Siler estimation failed for the following years: ",
    paste(failure_messages, collapse = "; ")
  )
}

fit_summary = do.call(
  rbind,
  lapply(fitted, function(item) {
    format_fit_summary(item$fit, item$year, sex)
  })
)

curve_rows = lapply(fitted, function(item) {
  year = item$year
  theta = item$fit$theta
  year_data = item$data
  components = siler_hazard_components(year_data$age, theta)
  data.frame(
    year = year,
    age = year_data$age,
    deaths = year_data$deaths,
    exposure = year_data$exposure,
    observed_mx = year_data$deaths / year_data$exposure,
    fitted_mx = components$juvenile + components$background + components$senescent,
    juvenile = components$juvenile,
    background = components$background,
    senescent = components$senescent,
    stringsAsFactors = FALSE
  )
})
fit_curves = do.call(rbind, curve_rows)

life_expectancy_rows = lapply(fitted, function(item) {
  year = item$year
  theta = item$fit$theta
  ages = if (year %in% analysis_config$benchmark_years) {
    analysis_config$life_expectancy_ages
  } else {
    analysis_config$table_ages
  }
  values = t(vapply(
    ages,
    function(age) {
      siler_remaining_life(
        age = age,
        theta = theta,
        absolute_tolerance = analysis_config$series_absolute_tolerance,
        relative_tolerance = analysis_config$quadrature_relative_tolerance
      )
    },
    numeric(7)
  ))
  values = as.data.frame(values)
  numerical_values = t(vapply(
    ages,
    function(age) {
      c(
        quadrature = siler_remaining_life_quadrature(
          age = age,
          theta = theta,
          relative_tolerance = analysis_config$quadrature_relative_tolerance
        ),
        siler_remaining_life_partial(
          age = age,
          theta = theta,
          maximum_indices = 0:3,
          relative_tolerance = analysis_config$quadrature_relative_tolerance
        )
      )
    },
    numeric(5)
  ))
  numerical_values = as.data.frame(numerical_values)
  data.frame(
    year = year,
    age = ages,
    siler_ex = values$exact,
    quadrature_ex = numerical_values$quadrature,
    partial_m0 = numerical_values$m0,
    partial_m1 = numerical_values$m1,
    partial_m2 = numerical_values$m2,
    partial_m3 = numerical_values$m3,
    gm_lower = values$gm_lower,
    gm_upper = values$gm_upper,
    truncation_bound = values$truncation_bound,
    series_terms = as.integer(values$terms),
    residual_juvenile_exposure = values$A,
    senescent_scale = values$G,
    stringsAsFactors = FALSE
  )
})
life_expectancy = do.call(rbind, life_expectancy_rows)
life_expectancy = merge(
  life_expectancy,
  life_table,
  by = c("year", "age"),
  all.x = TRUE,
  sort = TRUE
)
life_expectancy$model_minus_hmd = life_expectancy$siler_ex - life_expectancy$hmd_ex
life_expectancy$gm_relative_error_percent =
  100 * (life_expectancy$gm_upper - life_expectancy$siler_ex) /
  life_expectancy$siler_ex

quadrature_check = life_expectancy[
  life_expectancy$year %in% analysis_config$benchmark_years &
    life_expectancy$age %in% analysis_config$table_ages,
  c("year", "age", "siler_ex", "quadrature_ex", "truncation_bound"),
  drop = FALSE
]
names(quadrature_check)[names(quadrature_check) == "siler_ex"] = "series"
names(quadrature_check)[names(quadrature_check) == "quadrature_ex"] = "quadrature"
names(quadrature_check)[names(quadrature_check) == "truncation_bound"] =
  "analytic_error_bound"
quadrature_check$absolute_difference = abs(
  quadrature_check$series - quadrature_check$quadrature
)
quadrature_check = quadrature_check[
  c(
    "year", "age", "series", "quadrature", "absolute_difference",
    "analytic_error_bound"
  )
]

positive_observed = fit_curves$observed_mx > 0
residual_summary = aggregate(
  abs(log(fit_curves$observed_mx[positive_observed]) -
        log(fit_curves$fitted_mx[positive_observed])),
  by = list(year = fit_curves$year[positive_observed]),
  FUN = max
)
names(residual_summary)[2] = "maximum_absolute_log_rate_residual"
fit_summary = merge(fit_summary, residual_summary, by = "year", all.x = TRUE)

derived_directory = project_path("data", "derived")
dir.create(derived_directory, recursive = TRUE, showWarnings = FALSE)
results = list(
  config = analysis_config,
  source_files = c(
    deaths = death_file,
    exposures = exposure_file,
    life_table = life_table_file
  ),
  fit_objects = fitted,
  fit_summary = fit_summary,
  fit_curves = fit_curves,
  life_expectancy = life_expectancy,
  quadrature_check = quadrature_check
)

saveRDS(results, file.path(derived_directory, "sweden_siler_results.rds"))
write.csv(
  fit_summary,
  file.path(derived_directory, "sweden_siler_fit_summary.csv"),
  row.names = FALSE
)
write.csv(
  life_expectancy,
  file.path(derived_directory, "sweden_siler_life_expectancy.csv"),
  row.names = FALSE
)
write.csv(
  quadrature_check,
  file.path(derived_directory, "sweden_siler_quadrature_check.csv"),
  row.names = FALSE
)

cat("Saved fitted results to data/derived/sweden_siler_results.rds\n")
cat("Review optimizer convergence codes before generating outputs.\n")
print(fit_summary)
print(quadrature_check)

