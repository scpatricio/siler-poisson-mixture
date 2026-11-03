# Run the complete empirical workflow from any working directory.

r_directory = file.path(getwd(), "R")
source(file.path(r_directory, "config.R"))

required_packages = c("expint", "ggplot2")
missing_packages = required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install the required packages before running the analysis: ",
    paste(missing_packages, collapse = ", ")
  )
}

required_files = c(
  file.path(r_directory, "config.R"),
  file.path(r_directory, "siler_functions.R"),
  project_path("data", "hmd", "Deaths_1x1.txt"),
  project_path("data", "hmd", "Exposures_1x1.txt"),
  project_path("data", "hmd", "fltper_1x1.txt")
)
missing_files = required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop(
    "Required project files are missing: ",
    paste(missing_files, collapse = ", ")
  )
}

source(file.path(r_directory, "01_fit_sweden.R"))
source(file.path(r_directory, "02_make_figures.R"))
source(file.path(r_directory, "03_make_appendix_table.R"))
source(file.path(r_directory, "04_make_appendix_error_figure.R"))
source(file.path(r_directory, "05_make_immaturity_figure.R"))

cat("Analysis completed.\n")
cat("Review data/derived/sweden_siler_fit_summary.csv.\n")
cat("Figures: figures/\n")
cat("Appendix table: tables/table_numerical_verification.tex\n")
cat("Appendix figure: figures/figure_maximum_approximation_difference.pdf\n")

