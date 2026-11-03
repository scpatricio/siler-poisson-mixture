# Generate the appendix table that verifies the exact series numerically.

r_directory = file.path(getwd(), "R")
source(file.path(r_directory, "config.R"))

results_path = project_path("data", "derived", "sweden_siler_results.rds")
if (!file.exists(results_path)) {
  stop("Run R/01_fit_sweden.R before generating the appendix table.")
}

quadrature_check = readRDS(results_path)$quadrature_check
quadrature_check = quadrature_check[
  order(quadrature_check$year, quadrature_check$age),
  ,
  drop = FALSE
]

table_directory = project_path("tables")
dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)
table_path = file.path(table_directory, "table_numerical_verification.tex")

format_years = function(value) sprintf("%.4f", value)
format_scientific = function(value) {
  if (!is.finite(value)) {
    return("--")
  }
  if (abs(value) < 1e-16) {
    return("$<10^{-16}$")
  }
  exponent = floor(log10(abs(value)))
  coefficient = value / 10^exponent
  sprintf("$%.2f\\times 10^{%d}$", coefficient, exponent)
}

rows = vapply(
  seq_len(nrow(quadrature_check)),
  function(index) {
    row = quadrature_check[index, ]
    paste(
      row$year,
      row$age,
      format_years(row$series),
      format_scientific(row$absolute_difference),
      format_scientific(row$analytic_error_bound),
      sep = " & "
    )
  },
  character(1)
)
rows = paste0(rows, " \\\\")

latex = c(
  "\\begin{table}[htb!]",
  "\\centering",
  "\\small",
  paste0(
    "\\caption{Numerical verification of the exact Siler series. ",
    "The value column reports the series to four decimal places; adaptive ",
    "quadrature gives the same displayed values. The absolute difference ",
    "compares the series, stopped at an absolute truncation tolerance of ",
    "$10^{-12}$ years, with adaptive quadrature. The final column is the ",
    "analytic bound on the omitted series tail. Values below $10^{-16}$ ",
    "are printed as $<10^{-16}$. Differences can exceed the bound because ",
    "they also include the numerical error of the special-function routine ",
    "and of the quadrature, which the bound does not cover.}"
  ),
  "\\label{tab:numerical-verification}",
  "\\begin{tabular}{rrrrr}",
  "\\toprule",
  "Year & Age & Value & Absolute difference & Analytic bound \\\\",
  "\\midrule",
  rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)

writeLines(latex, table_path, useBytes = TRUE)
cat("Saved appendix table to tables/table_numerical_verification.tex\n")
