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
  c("year", "siler_ex", "gm_upper", "remaining_immaturity_hazard")
]
birth = merge(birth, fit_summary[c("year", "a2", "b1")], by = "year")
birth = birth[order(birth$year), ]

# A_0 is the cumulative immaturity hazard at birth; gm_upper is the Gompertz--Makeham
# life expectancy with the immaturity hazard removed.
birth$immaturity_probability = 1 - exp(-birth$remaining_immaturity_hazard)
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

# Numbers reported in Section 5.1: the loss at birth in selected years, the
# cause-deleted decomposition of the gain in e0 between the first and last
# fitted years, and the accuracy of the leading-term approximation.
first_year = min(birth$year)
last_year = max(birth$year)
first = birth[birth$year == first_year, ]
last = birth[birth$year == last_year, ]

gain_e0 = last$siler_ex - first$siler_ex
gain_without_immaturity = last$gm_upper - first$gm_upper
decline_in_loss = first$loss - last$loss
share_of_gain = decline_in_loss / gain_e0

loss_ratio = birth$loss / birth$loss_leading_term
leading_term_gap = birth$loss_leading_term - birth$loss

decomposition = data.frame(
  first_year = first_year,
  last_year = last_year,
  loss_first_year = first$loss,
  loss_1950 = if (1950 %in% birth$year) birth$loss[birth$year == 1950] else NA,
  loss_last_year = last$loss,
  gain_e0 = gain_e0,
  gain_without_immaturity = gain_without_immaturity,
  decline_in_loss = decline_in_loss,
  share_of_gain = share_of_gain,
  minimum_loss_ratio = min(loss_ratio),
  maximum_loss_ratio = max(loss_ratio),
  maximum_leading_term_gap = max(leading_term_gap)
)

write.csv(
  decomposition,
  project_path("data", "derived", "sweden_siler_immaturity_decomposition.csv"),
  row.names = FALSE
)

cat(sprintf(
  "Loss at birth: %.2f years in %d, %.2f in 1950, %.3f in %d.\n",
  first$loss, first_year, decomposition$loss_1950, last$loss, last_year
))
cat(sprintf(
  "Gain in e0: %.1f years; gain with the immaturity hazard removed: %.1f years.\n",
  gain_e0, gain_without_immaturity
))
cat(sprintf(
  "Decline of the loss: %.2f years, or %.1f%% of the gain in e0.\n",
  decline_in_loss, 100 * share_of_gain
))
cat(sprintf(
  "Loss as a share of its upper bound: %.1f%% to %.2f%%; maximum gap %.3f years.\n",
  100 * min(loss_ratio), 100 * max(loss_ratio), max(leading_term_gap)
))

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
plot_theme = ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(color = "#E5E5E5", linewidth = 0.3),
    strip.text = ggplot2::element_text(face = "italic", size = 11),
    panel.spacing = grid::unit(1.1, "lines"),
    axis.title = ggplot2::element_text(face = "plain"),
    legend.position = "bottom",
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 10.6),
    plot.title.position = "plot",
    plot.title = ggplot2::element_text(face = "bold", size = 15),
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
      "Exact loss" = colors[["immaturity"]],
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
    x = "Year",
    y = NULL
  ) +
  plot_theme

figure_directory = project_path("figures")
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(
  filename = file.path(figure_directory, "figure_immaturity_loss.pdf"),
  plot = figure_immaturity,
  width = 9,
  height = 5,
  units = "in",
  device = grDevices::cairo_pdf
)

cat("Saved figure_immaturity_loss.pdf, sweden_siler_immaturity_loss.csv, and sweden_siler_immaturity_decomposition.csv.\n")
