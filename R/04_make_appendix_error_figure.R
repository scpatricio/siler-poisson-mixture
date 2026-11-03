# Plot worst-case fixed-m approximation differences across all fitted years.

r_directory = file.path(getwd(), "R")
source(file.path(r_directory, "config.R"))
source(file.path(r_directory, "siler_functions.R"))

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required. Install it with install.packages('ggplot2').")
}

results_path = project_path("data", "derived", "sweden_siler_results.rds")
if (!file.exists(results_path)) {
  stop("Run R/01_fit_sweden.R before generating the appendix figure.")
}
results = readRDS(results_path)
fit_summary = results$fit_summary[
  results$fit_summary$year %in% analysis_config$years,
  ,
  drop = FALSE
]

appendix_ages = 0:100
maximum_indices = 0:4

calculate_year_differences = function(year) {
  fitted = fit_summary[fit_summary$year == year, , drop = FALSE]
  if (nrow(fitted) != 1L) {
    stop("Expected one fitted parameter row for year ", year, ".")
  }
  theta = unlist(fitted[1L, c("a1", "b1", "a2", "a3", "b3")])
  do.call(
    rbind,
    lapply(appendix_ages, function(age) {
      fitted_life = siler_remaining_life_quadrature(
        age = age,
        theta = theta,
        relative_tolerance = analysis_config$quadrature_relative_tolerance
      )
      absolute_difference = siler_remaining_life_truncation_difference(
        age = age,
        theta = theta,
        maximum_indices = maximum_indices,
        tail_tolerance = analysis_config$poisson_tail_tolerance,
        relative_tolerance = analysis_config$quadrature_relative_tolerance
      )
      data.frame(
        year = year,
        age = age,
        m = maximum_indices,
        absolute_difference = unname(absolute_difference),
        absolute_relative_difference_percent =
          100 * unname(absolute_difference) / fitted_life
      )
    })
  )
}

years = analysis_config$years
number_of_cores = min(length(years), analysis_config$parallel_cores)
if (.Platform$OS.type == "unix" && number_of_cores > 1L) {
  yearly_differences = parallel::mclapply(
    years,
    calculate_year_differences,
    mc.cores = number_of_cores,
    mc.preschedule = FALSE
  )
} else {
  yearly_differences = lapply(years, calculate_year_differences)
}
failed = vapply(yearly_differences, inherits, logical(1), "try-error")
if (any(failed)) {
  stop(
    "Approximation differences failed for years: ",
    paste(years[failed], collapse = ", ")
  )
}
difference_data = do.call(rbind, yearly_differences)

maximum_absolute = aggregate(
  absolute_difference ~ age + m,
  data = difference_data,
  FUN = max
)
maximum_relative = aggregate(
  absolute_relative_difference_percent ~ age + m,
  data = difference_data,
  FUN = max
)
maximum_differences = merge(
  maximum_absolute,
  maximum_relative,
  by = c("age", "m"),
  sort = TRUE
)

derived_directory = project_path("data", "derived")
dir.create(derived_directory, recursive = TRUE, showWarnings = FALSE)
write.csv(
  maximum_differences,
  file.path(derived_directory, "sweden_siler_maximum_approximation_differences.csv"),
  row.names = FALSE
)

plot_data = rbind(
  data.frame(
    age = maximum_differences$age,
    m = maximum_differences$m,
    measure = "Maximum absolute difference (years)",
    value = maximum_differences$absolute_difference
  ),
  data.frame(
    age = maximum_differences$age,
    m = maximum_differences$m,
    measure = "Maximum absolute relative difference (%)",
    value = maximum_differences$absolute_relative_difference_percent
  )
)
plot_data = plot_data[plot_data$value >= 1e-50, , drop = FALSE]
plot_data$series = factor(
  paste0("m = ", plot_data$m),
  levels = paste0("m = ", maximum_indices)
)
plot_data$measure = factor(
  plot_data$measure,
  levels = c(
    "Maximum absolute difference (years)",
    "Maximum absolute relative difference (%)"
  )
)

appendix_figure = ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = age,
    y = value,
    color = series,
    group = series
  )
) +
  ggplot2::geom_line(linewidth = 1.0, na.rm = TRUE) +
  ggplot2::facet_wrap(~measure, nrow = 1L, scales = "free_y") +
  ggplot2::scale_color_manual(
    values = c(
      "m = 0" = "#D55E00",
      "m = 1" = "#E69F00",
      "m = 2" = "#009E73",
      "m = 3" = "#0072B2",
      "m = 4" = "#CC79A7"
    )
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq(0, 100, by = 20),
    limits = c(0, 100),
    expand = ggplot2::expansion(mult = c(0.01, 0.01))
  ) +
  ggplot2::scale_y_log10(
    limits = c(1e-50, 1),
    breaks = 10^c(0, -10, -20, -30, -40, -50),
    labels = scales::label_log(base = 10)
  ) +
  ggplot2::labs(
    title = "Maximum fixed-m approximation differences by age",
    subtitle = expression(
      paste(
        "Maximum at each age across Swedish female mortality schedules, ",
        "1900-2020; values below ", 10^-50, " are not displayed"
      )
    ),
    x = "Age",
    y = NULL
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(nrow = 1L, byrow = TRUE)
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(color = "#E5E5E5", linewidth = 0.3),
    strip.text = ggplot2::element_text(face = "bold"),
    panel.spacing = grid::unit(1.1, "lines"),
    legend.position = "bottom",
    legend.title = ggplot2::element_blank(),
    plot.title.position = "plot",
    plot.title = ggplot2::element_text(face = "bold", size = 14),
    plot.subtitle = ggplot2::element_text(size = 10.5, color = "#3F3F3F"),
    plot.margin = ggplot2::margin(10, 12, 8, 10)
  )

figure_directory = project_path("figures")
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(
  filename = file.path(
    figure_directory,
    "figure_maximum_approximation_difference.pdf"
  ),
  plot = appendix_figure,
  width = 10.5,
  height = 5.2,
  units = "in",
  device = grDevices::cairo_pdf
)

cat("Saved appendix maximum-difference figure to figures/.\n")
