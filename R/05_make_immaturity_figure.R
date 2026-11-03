# Life expectancy at birth lost to immaturity mortality in the Swedish Siler
# fits, with its leading-term decomposition. Uses the derived results from
# 01_fit_sweden.R and writes one figure and one CSV file.

r_directory = file.path(getwd(), "R")
source(file.path(r_directory, "config.R"))

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required. Install it with install.packages('ggplot2').")
}

life_expectancy_path = project_path(
  "data",
  "derived",
  "sweden_siler_life_expectancy.csv"
)
fit_summary_path = project_path("data", "derived", "sweden_siler_fit_summary.csv")
if (!file.exists(life_expectancy_path) || !file.exists(fit_summary_path)) {
  stop("Run R/01_fit_sweden.R before generating the immaturity figure.")
}

life_expectancy = read.csv(life_expectancy_path)
fit_summary = read.csv(fit_summary_path)

birth = life_expectancy[
  life_expectancy$age == 0 & life_expectancy$year %in% analysis_config$years,
  c("year", "siler_ex", "gm_upper", "residual_juvenile_exposure")
]
birth = merge(birth, fit_summary[c("year", "a2", "b1")], by = "year")
birth = birth[order(birth$year), ]

# A_0 is the total immaturity hazard; gm_upper is the Gompertz--Makeham
# life expectancy with the immaturity hazard removed.
birth$immaturity_probability = 1 - exp(-birth$residual_juvenile_exposure)
birth$loss = birth$gm_upper - birth$siler_ex
birth$loss_leading_term = birth$immaturity_probability * birth$gm_upper
birth$loss_lower_bound = birth$immaturity_probability *
  (birth$gm_upper - 1 / (birth$a2 + birth$b1))

if (any(birth$loss > birth$loss_leading_term + 1e-9) ||
    any(birth$loss < birth$loss_lower_bound - 1e-9)) {
  stop("The immaturity loss violates its analytic bounds.")
}

write.csv(
  birth,
  project_path("data", "derived", "sweden_siler_immaturity_loss.csv"),
  row.names = FALSE
)

panel_levels = c(
  "Life expectancy at birth (years)",
  "Life expectancy lost to immaturity mortality (years)"
)
series_levels = c(
  "Siler",
  "Immaturity hazard removed",
  "Exact loss",
  "Leading term"
)
plot_data = rbind(
  data.frame(
    year = birth$year,
    value = birth$siler_ex,
    series = "Siler",
    panel = panel_levels[1]
  ),
  data.frame(
    year = birth$year,
    value = birth$gm_upper,
    series = "Immaturity hazard removed",
    panel = panel_levels[1]
  ),
  data.frame(
    year = birth$year,
    value = birth$loss,
    series = "Exact loss",
    panel = panel_levels[2]
  ),
  data.frame(
    year = birth$year,
    value = birth$loss_leading_term,
    series = "Leading term",
    panel = panel_levels[2]
  )
)
plot_data$panel = factor(plot_data$panel, levels = panel_levels)
plot_data$series = factor(plot_data$series, levels = series_levels)

colors = analysis_config$colors
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

figure_immaturity = ggplot2::ggplot(
  plot_data,
  ggplot2::aes(x = year, y = value, color = series, linetype = series)
) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::facet_wrap(~panel, ncol = 2, scales = "free_y") +
  ggplot2::scale_color_manual(
    values = c(
      "Siler" = colors[["fitted"]],
      "Immaturity hazard removed" = colors[["senescent"]],
      "Exact loss" = colors[["juvenile"]],
      "Leading term" = "#111111"
    )
  ) +
  ggplot2::scale_linetype_manual(
    values = c(
      "Siler" = "solid",
      "Immaturity hazard removed" = "solid",
      "Exact loss" = "solid",
      "Leading term" = "dashed"
    )
  ) +
  ggplot2::labs(
    title = "Life expectancy lost to immaturity mortality",
    subtitle = "Swedish females, 1900-2020; fitted Siler schedules",
    x = "Year",
    y = NULL
  ) +
  plot_theme

figure_directory = project_path("figures")
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(
  filename = file.path(figure_directory, "figure_immaturity_loss.pdf"),
  plot = figure_immaturity,
  width = 10.5,
  height = 4.6,
  units = "in",
  device = grDevices::cairo_pdf
)

cat("Saved figure_immaturity_loss.pdf and sweden_siler_immaturity_loss.csv.\n")
