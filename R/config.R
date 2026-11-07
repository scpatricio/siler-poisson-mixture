# Analysis configuration for the Swedish Siler application.

r_directory = file.path(getwd(), "R")
if (!file.exists(file.path(r_directory, "config.R"))) {
  stop(
    "Set the working directory to the repository root before running the ",
    "scripts, for example by opening siler-poisson-mixture.Rproj."
  )
}
project_directory = dirname(r_directory)
project_path = function(...) file.path(project_directory, ...)

siler_detected_cores = parallel::detectCores(logical = FALSE)
if (!is.finite(siler_detected_cores)) {
  siler_detected_cores = 1L
}

analysis_config = list(
  sex = "Female",
  years = 1900:2020,
  benchmark_years = seq(1920L, 2020L, by = 20L),
  minimum_age = 0L,
  maximum_age = 109L,
  table_ages = c(0L, 30L, 60L),
  life_expectancy_ages = c(0L, 30L, 60L),
  number_of_starts = 16L,
  optimizer_max_iterations = 5000L,
  series_absolute_tolerance = 1e-12,
  poisson_tail_tolerance = 1e-14,
  quadrature_relative_tolerance = 1e-10,
  parallel_cores = max(1L, siler_detected_cores - 1L),
  random_seed = 20261103L,
  colors = c(
    observed = "#6B6B6B",
    fitted = "#0072B2",
    immaturity = "#E69F00",
    background = "#009E73",
    senescent = "#CC79A7",
    reference = "#D55E00"
  )
)

