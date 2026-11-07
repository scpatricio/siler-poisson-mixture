# Siler remaining life expectancy as a Poisson mixture

Code and data to reproduce the results in

> Patricio, S. C. (2026). An exact Poisson mixture representation of Siler
> remaining life expectancy. *arXiv:2611.06376*.
> <https://doi.org/10.48550/arXiv.2611.06376>

```bibtex
@article{patricio2026siler,
  title   = {An Exact Poisson Mixture Representation of {Siler} Remaining Life Expectancy},
  author  = {Patricio, Silvio C.},
  journal = {arXiv preprint arXiv:2611.06376},
  year    = {2026},
  doi     = {10.48550/arXiv.2611.06376}
}
```

The paper shows that the remaining-lifetime distribution under the Siler
mortality model is Gompertz–Makeham with a Poisson-distributed Makeham term.
A latent-clock construction explains the representation: individuals carry a
Poisson number of constant-rate clocks, and selection removes those carrying
more clocks first. Integrating the mixture gives an exact series of
generalized exponential integrals for remaining life expectancy with a
closed-form truncation bound, and a two-sided bound on the life expectancy
lost to immaturity mortality. The application fits the Siler model to Swedish
female mortality, 1900–2020.

## Requirements

- R (≥ 4.1; the scripts use the native pipe `|>`)
- R packages `expint` and `ggplot2` (`scales` is installed with `ggplot2`)

```r
install.packages(c("expint", "ggplot2"))
```

## Reproducing the results

Set the working directory to the repository root (opening
`siler-poisson-mixture.Rproj` in RStudio does this) and run

```r
source("R/run_all.R")
```

or, from a terminal in the repository root,

```sh
Rscript R/run_all.R
```

The workflow fits 121 years and evaluates the series at ages 0–100 in every
year. On Unix-like systems the years are processed in parallel; the number of
cores is set in `R/config.R`.

## Repository structure

```
R/
  config.R                         analysis settings (years, ages, tolerances, cores)
  siler_functions.R                likelihood, optimizer, series and quadrature
  01_fit_sweden.R                  Poisson maximum-likelihood fits, 1900–2020
  02_make_figures.R                remaining life expectancy and series-accuracy figures
  03_make_appendix_table.R         numerical-verification table
  04_make_appendix_error_figure.R  maximum truncation error over all years
  05_make_immaturity_figure.R      life expectancy lost to immaturity mortality
  run_all.R                        runs scripts 01–05 in order
  validate_series.R                standalone check of the series on test parameters
data/
  hmd/                             Human Mortality Database input files (see below)
  derived/                         fitted parameters and intermediate results
figures/                           figures, as PDF
tables/                            LaTeX table
```

## Outputs

| File | Content |
|---|---|
| `figures/figure_remaining_life.pdf` | Fitted Siler and HMD remaining life expectancy at ages 0, 30 and 60 |
| `figures/figure_immaturity_loss.pdf` | Life expectancy at birth lost to immaturity mortality and its leading-term value |
| `figures/figure_series_convergence.pdf` | Truncation error of the series by age, 1950 |
| `figures/figure_maximum_approximation_difference.pdf` | Largest truncation error by age over all years |
| `tables/table_numerical_verification.tex` | Series versus adaptive quadrature in six benchmark years |
| `data/derived/sweden_siler_fit_summary.csv` | Estimated Siler parameters and convergence diagnostics by year |
| `data/derived/sweden_siler_life_expectancy.csv` | Series, partial sums, quadrature and HMD values at ages 0, 30 and 60 |
| `data/derived/sweden_siler_immaturity_loss.csv` | Loss to immaturity mortality at birth and its bounds, by year |
| `data/derived/sweden_siler_immaturity_decomposition.csv` | Numbers reported in Section 5.1: loss at birth, cause-deleted decomposition of the gain in e0, and accuracy of the leading term |
| `data/derived/sweden_siler_maximum_approximation_differences.csv` | Maximum truncation error by age and number of terms |
| `data/derived/sweden_siler_quadrature_check.csv` | Values reported in the verification table |
| `data/derived/sweden_siler_results.rds` | All fitted results, read by scripts 02–04 |

The derived files and the figures in this repository are the ones reported
in the paper. Running `R/run_all.R` regenerates them; the full run takes about
a minute on a laptop. Parameter estimates, life expectancies, truncation bounds
and approximation errors reproduce exactly. The "Absolute difference" column
of the verification table compares two numerical routines at the level of
floating-point error (about 1e-14 to 1e-11 years), so its last
digits can differ across platforms and library versions.

## Data

`data/hmd/` contains three files for Swedish females from the Human Mortality
Database: `Deaths_1x1.txt`, `Exposures_1x1.txt` and `fltper_1x1.txt` (HMD
version last modified 29 May 2025, Methods Protocol v6). The files are
unchanged. Data constructed by the HMD are distributed under the CC BY 4.0
license. Please cite them as

> HMD. Human Mortality Database. Max Planck Institute for Demographic Research
> (Germany), University of California, Berkeley (USA), and French Institute
> for Demographic Studies (France). Available at www.mortality.org.

To use a newer release, download the same three files for Sweden from
<https://www.mortality.org/> and replace those in `data/hmd/`.

## License

The code is released under the MIT License (see `LICENSE`). The HMD data
remain under the CC BY 4.0 license of the Human Mortality Database.
