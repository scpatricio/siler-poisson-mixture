# Generate publication-quality figures from the fitted Swedish Siler results.
# Outputs are written as vector PDF files only.

r_directory = file.path(getwd(), "R")
source(file.path(r_directory, "config.R"))
source(file.path(r_directory, "siler_functions.R"))

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required. Install it with install.packages('ggplot2').")
}

results_path = project_path("data", "derived", "sweden_siler_results.rds")
if (!file.exists(results_path)) {
  stop("Run R/01_fit_sweden.R before generating figures.")
}
results = readRDS(results_path)
fit_summary = results$fit_summary
fit_curves = results$fit_curves
life_expectancy = results$life_expectancy
fit_summary = fit_summary[
  fit_summary$year %in% analysis_config$years,
  ,
  drop = FALSE
]
fit_curves = fit_curves[
  fit_curves$year %in% analysis_config$years,
  ,
  drop = FALSE
]
life_expectancy = life_expectancy[
  life_expectancy$year %in% analysis_config$years,
  ,
  drop = FALSE
]
colors = analysis_config$colors

if (any(fit_summary$convergence != 0)) {
  warning(
    "At least one optimizer did not report convergence. Review fit diagnostics ",
    "before using the figures in the manuscript."
  )
}
figure_directory = project_path("figures")
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)

year_levels = as.character(analysis_config$years)
fit_curves$year_label = factor(fit_curves$year, levels = analysis_config$years)
life_expectancy$year_label = factor(
  life_expectancy$year,
  levels = analysis_config$years
)

benchmark_curves = fit_curves[
  fit_curves$year %in% analysis_config$benchmark_years,
  ,
  drop = FALSE
]
benchmark_curves$year_label = factor(
  benchmark_curves$year,
  levels = analysis_config$benchmark_years
)
plot_theme = ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(color = "#E5E5E5", linewidth = 0.3),
    strip.text = ggplot2::element_text(face = "bold"),
    panel.spacing = grid::unit(1.1, "lines"),
    axis.title = ggplot2::element_text(face = "plain"),
    legend.position = "bottom",
    legend.title = ggplot2::element_blank(),
    plot.title.position = "plot",
    plot.title = ggplot2::element_text(face = "bold", size = 14),
    plot.subtitle = ggplot2::element_text(size = 10.5, color = "#3F3F3F"),
    plot.margin = ggplot2::margin(10, 12, 8, 10)
  )

observed_positive = benchmark_curves[
  benchmark_curves$observed_mx > 0,
  ,
  drop = FALSE
]

figure_mortality_fit = ggplot2::ggplot() +
  ggplot2::geom_point(
    data = observed_positive,
    ggplot2::aes(x = age, y = observed_mx, color = "Observed"),
    # size = 1,
    colour = 'gray60',
    alpha = .5
  ) +
  ggplot2::geom_line(
    data = benchmark_curves,
    ggplot2::aes(x = age, y = fitted_mx, color = "Fitted Siler"),
    linewidth = .8
  ) +
  ggplot2::facet_wrap(~year_label, ncol = 3) +
  ggplot2::scale_color_manual(
    values = c(
      "Observed" = unname(colors["observed"]),
      "Fitted Siler" = unname(colors["fitted"])
    ),
    breaks = c("Observed", "Fitted Siler")
  ) +
  ggplot2::scale_y_log10(
    limits = c(1e-5, 10),
    breaks = 10^(-5:1),
    labels = scales::label_log(base = 10)
  ) +
  ggplot2::scale_x_continuous(
    breaks = c(0, 20, 40, 60, 80, 100),
    limits = c(0, analysis_config$maximum_age),
    expand = ggplot2::expansion(mult = c(0.01, 0.02))
  ) +
  ggplot2::labs(
    title = "Siler mortality fits for Swedish females",
    subtitle = "Observed death rates and Poisson maximum-likelihood fits",
    x = "Age (years)",
    y = "Mortality rate, log scale",
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(order = 1, nrow = 1, byrow = TRUE)
  ) +
  plot_theme

selected_life_plot = life_expectancy[
  life_expectancy$age %in% analysis_config$table_ages,
  ,
  drop = FALSE
]

required_life_columns = c(
  "hmd_ex", "quadrature_ex"
)

if (!all(required_life_columns %in% names(selected_life_plot))) {
  stop("Run the revised R/01_fit_sweden.R before generating figures.")
}

if (any(!is.finite(selected_life_plot$hmd_ex))) {
  stop("HMD remaining life expectancy is missing for at least one plotted year-age pair.")
}

remaining_value_data = rbind(
  data.frame(
    year = selected_life_plot$year,
    age = selected_life_plot$age,
    series = "HMD",
    value = selected_life_plot$hmd_ex
  ),
  data.frame(
    year = selected_life_plot$year,
    age = selected_life_plot$age,
    series = "Fitted Siler",
    value = selected_life_plot$quadrature_ex
  )
)

remaining_value_data$measure = "LE (years)"

fitted_minus_hmd = selected_life_plot$quadrature_ex - selected_life_plot$hmd_ex

remaining_error_data = rbind(
  data.frame(
    year = selected_life_plot$year,
    age = selected_life_plot$age,
    measure = "Difference (years)",
    series = "Fitted Siler",
    value = fitted_minus_hmd
  ),
  data.frame(
    year = selected_life_plot$year,
    age = selected_life_plot$age,
    measure = "Relative difference (%)",
    series = "Fitted Siler",
    value = 100 * fitted_minus_hmd / selected_life_plot$hmd_ex
  )
)

remaining_plot_data = rbind(remaining_value_data, remaining_error_data)

remaining_plot_data$age_label = paste("Age", remaining_plot_data$age)

measure_levels = c(
  "LE (years)",
  "Difference (years)",
  "Relative difference (%)"
)

remaining_plot_data$measure = factor(
  remaining_plot_data$measure,
  levels = measure_levels
)

remaining_plot_data$age_label = factor(
  remaining_plot_data$age_label,
  levels = paste("Age", analysis_config$table_ages)
)

reference_rows = expand.grid(
  age_label = factor(
    paste("Age", analysis_config$table_ages),
    levels = paste("Age", analysis_config$table_ages)
  ),
  measure = factor(measure_levels[-1], levels = measure_levels),
  KEEP.OUT.ATTRS = FALSE
)

reference_rows$yintercept = 0

figure_remaining_life = ggplot2::ggplot(
  remaining_plot_data,
  ggplot2::aes(x = year, y = value, color = series)
) +
  ggplot2::geom_hline(
    data = reference_rows,
    ggplot2::aes(yintercept = yintercept),
    inherit.aes = FALSE,
    color = "#777777",
    linewidth = 0.4
  ) +
  ggplot2::geom_line(ggplot2::aes(group = series), linewidth = 0.8) +
  ggplot2::geom_line(
    data = remaining_plot_data[
      remaining_plot_data$series == "Fitted Siler",
      ,
      drop = FALSE
    ],
    ggplot2::aes(group = series),
    color = "#0072B2",
    linetype = "dashed",
    linewidth = 0.95,
    show.legend = FALSE
  ) +
  ggplot2::facet_grid(measure ~ age_label, scales = "free_y") +
  ggplot2::scale_color_manual(
    values = c(
      "HMD" = "#111111",
      "Fitted Siler" = "#0072B2"
    ),
    breaks = c("HMD", "Fitted Siler")
  ) +
  ggplot2::scale_x_continuous(breaks = seq(1900, 2020, by = 20)) +
  ggplot2::labs(
    title = "Siler remaining life expectancy in Sweden, 1900-2020",
    subtitle = "HMD values and fitted Siler values for Swedish females",
    x = "Year",
    y = NULL
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(nrow = 1, byrow = TRUE)
  ) +
  plot_theme

performance_year = 1950L

performance_ages = 0:100

performance_fit = fit_summary[
  fit_summary$year == performance_year,
  ,
  drop = FALSE
]

if (nrow(performance_fit) != 1L) {
  stop(
    "The approximation illustration requires exactly one fitted row for ",
    performance_year,
    "."
  )
}

performance_theta = unlist(
  performance_fit[1L, c("a1", "b1", "a2", "a3", "b3")],
  use.names = TRUE
)

performance_quadrature = vapply(
  performance_ages,
  function(age) {
    siler_remaining_life_quadrature(
      age = age,
      theta = performance_theta,
      relative_tolerance = analysis_config$quadrature_relative_tolerance
    )
  },
  numeric(1)
)

performance_absolute_differences = t(vapply(
  performance_ages,
  function(age) {
    siler_remaining_life_truncation_difference(
      age = age,
      theta = performance_theta,
      maximum_indices = 0:4,
      tail_tolerance = analysis_config$poisson_tail_tolerance,
      relative_tolerance = analysis_config$quadrature_relative_tolerance
    )
  },
  numeric(5)
))

colnames(performance_absolute_differences) = paste0("m", 0:4)

performance_error_data = do.call(
  rbind,
  lapply(0:4, function(index) {
    absolute_difference =
      performance_absolute_differences[, paste0("m", index)]
    displayed_absolute_difference = ifelse(
      absolute_difference >= 1e-50,
      absolute_difference,
      NA_real_
    )
    rbind(
      data.frame(
        age = performance_ages,
        series = paste0("m = ", index),
        measure = "Absolute difference (years)",
        value = displayed_absolute_difference
      ),
      data.frame(
        age = performance_ages,
        series = paste0("m = ", index),
        measure = "Absolute relative difference (%)",
        value = 100 * displayed_absolute_difference / performance_quadrature
      )
    )
  })
)

performance_error_data = performance_error_data[
  performance_error_data$age <= 50,
  ,
  drop = FALSE
]

approximation_colors = c(
  "m = 0" = "#D55E00",
  "m = 1" = "#E69F00",
  "m = 2" = "#009E73",
  "m = 3" = "#0072B2",
  "m = 4" = "#CC79A7"
)

approximation_breaks = paste0("m = ", 0:4)

age_scale = ggplot2::scale_x_continuous(
  breaks = seq(0, 50, by = 10),
  limits = c(0, 50),
  expand = ggplot2::expansion(mult = c(0.01, 0.01))
)

approximation_panel_theme = plot_theme +
  ggplot2::theme(
    legend.position = "none",
    plot.title = ggplot2::element_text(face = "bold", size = 10.5, hjust = 0.5),
    plot.margin = ggplot2::margin(3, 5, 3, 5)
  )

absolute_difference_data = performance_error_data[
  performance_error_data$measure == "Absolute difference (years)",
  ,
  drop = FALSE
]

figure_age_absolute_difference = ggplot2::ggplot(
  absolute_difference_data,
  ggplot2::aes(x = age, y = value, color = series)
) +
  ggplot2::geom_line(
    ggplot2::aes(group = series),
    linewidth = 1.0,
    na.rm = TRUE
  ) +
  ggplot2::scale_color_manual(
    values = approximation_colors,
    breaks = approximation_breaks
  ) +
  ggplot2::scale_y_log10(
    limits = c(1e-50, 1),
    breaks = 10^c(0, -10, -20, -30, -40, -50),
    labels = scales::label_log(base = 10)
  ) +
  age_scale +
  ggplot2::labs(
    title = "Absolute difference",
    x = "Age",
    y = "Years, log scale"
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(nrow = 1L, byrow = TRUE)
  ) +
  approximation_panel_theme

absolute_relative_difference_data = performance_error_data[
  performance_error_data$measure == "Absolute relative difference (%)",
  ,
  drop = FALSE
]

figure_age_absolute_relative_difference = ggplot2::ggplot(
  absolute_relative_difference_data,
  ggplot2::aes(x = age, y = value, color = series)
) +
  ggplot2::geom_line(
    ggplot2::aes(group = series),
    linewidth = 1.0,
    na.rm = TRUE
  ) +
  ggplot2::scale_color_manual(
    values = approximation_colors,
    breaks = approximation_breaks
  ) +
  ggplot2::scale_y_log10(
    limits = c(1e-50, 1),
    breaks = 10^c(0, -10, -20, -30, -40, -50),
    labels = scales::label_log(base = 10)
  ) +
  age_scale +
  ggplot2::labs(
    title = "Absolute relative difference",
    x = "Age",
    y = "Percent, log scale"
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(nrow = 1L, byrow = TRUE)
  ) +
  approximation_panel_theme

extract_plot_legend = function(plot) {
  plot_table = ggplot2::ggplotGrob(
    plot + ggplot2::theme(
      legend.position = "bottom",
      legend.box = "horizontal"
    )
  )
  legend_indices = which(vapply(
    plot_table$grobs,
    function(grob) grepl("^guide-box", grob$name) &&
      !inherits(grob, "zeroGrob"),
    logical(1)
  ))
  if (length(legend_indices) == 0L) {
    stop("Could not extract the shared approximation legend.")
  }
  plot_table$grobs[[legend_indices[1L]]]
}

approximation_legend = extract_plot_legend(figure_age_absolute_difference)

draw_age_approximation_figure = function() {
  grid::grid.newpage()
  figure_layout = grid::grid.layout(
    nrow = 4L,
    ncol = 2L,
    heights = grid::unit(c(0.09, 0.07, 0.74, 0.10), "npc"),
    widths = grid::unit(rep(1, 2), "null")
  )
  grid::pushViewport(grid::viewport(layout = figure_layout))
  grid::grid.text(
    "Accuracy of the truncated Siler series by attained age",
    x = grid::unit(0.01, "npc"),
    just = "left",
    gp = grid::gpar(fontsize = 14, fontface = "bold"),
    vp = grid::viewport(layout.pos.row = 1L, layout.pos.col = 1:2)
  )
  grid::grid.text(
    substitute(
      paste(
        "Swedish females, ", year,
        "; omitted-tail errors below ", 10^-50,
        " are not displayed"
      ),
      list(year = performance_year)
    ),
    x = grid::unit(0.01, "npc"),
    just = "left",
    gp = grid::gpar(fontsize = 10.5, col = "#3F3F3F"),
    vp = grid::viewport(layout.pos.row = 2L, layout.pos.col = 1:2)
  )
  panel_plots = list(
    figure_age_absolute_difference,
    figure_age_absolute_relative_difference
  )
  for (column in seq_along(panel_plots)) {
    grid::pushViewport(
      grid::viewport(layout.pos.row = 3L, layout.pos.col = column)
    )
    grid::grid.draw(ggplot2::ggplotGrob(panel_plots[[column]]))
    grid::popViewport()
  }
  grid::pushViewport(
    grid::viewport(layout.pos.row = 4L, layout.pos.col = 1:2)
  )
  grid::grid.draw(approximation_legend)
  grid::popViewport()
  grid::popViewport()
}

save_figure = function(plot, stem, width = 7.2, height = 5.6) {
  ggplot2::ggsave(
    filename = file.path(figure_directory, paste0(stem, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf
  )
}

save_figure(figure_remaining_life, "figure_remaining_life", width = 10.5, height = 7.5)
grDevices::cairo_pdf(
  filename = file.path(figure_directory, "figure_series_convergence.pdf"),
  width = 10.5,
  height = 5.2
)
draw_age_approximation_figure()
grDevices::dev.off()

cat("Saved figures to figures/.\n")

