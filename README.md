<img src="man/figures/logo.png" width="30%" alt="recalc logo" />

# recalc

<!-- badges: start -->
[![R-CMD-check](https://github.com/ianhussey/recalc/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ianhussey/recalc/actions/workflows/R-CMD-check.yaml)

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21439896.svg)](https://doi.org/10.5281/zenodo.21439896)
<!-- badges: end -->

`recalc` is a forensic meta-science R package for conducting trustworthiness assessments on reported rounded summary statistics using recalculation checks. 

It recalculates *p*-values and standardised effect sizes (e.g. Cohen's *d*) from reported summary and test statistics, taking into account (a) the rounding of the inputs and outputs and (b) the analytic choices a paper usually leaves unstated (e.g. Student's vs Welch's *t*; whether reported SEs may have been reported as SDs).

## The idea

A reported statistical result is not a single number. It was computed from inputs that are rounded, and via an analysis whose choices a paper rarely states in full. So a reported result is compatible with a *range* of recalculated values, not one.

`recalc` recomputes the result across the multiverse of those choices — the rounding of every input and output, Student's vs Welch's *t*, one- vs two-sided tests, SD/SE confusion, and so on — and reports the resulting interval and whether the value printed in the paper falls inside it. A result that cannot be reproduced under *any* defensible choice is thereby distinguished from one that merely depends on an unstated convention.

As a recalculation method it is conceptually related to
[statcheck](https://statcheck.io/). It was inspired by:

- the `rpsychi` R package (v0.8, no longer on CRAN due to lack of maintenance),
- Nick Brown's `f_range()` function (<https://steamtraen.blogspot.com/2018/01/>).

## Installation

`recalc` is currently available from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("ianhussey/recalc")
```

## Usage

See the vignette (`vignette("recalc")`) for a worked introduction.

### Independent-samples *t*-test

Give the reported group means, SDs, and *N*s (with the number of decimals each was reported to). The function returns the range of *p*-values and Cohen's *d* values compatible with those inputs, and checks the reported *p* and *d* against it.

``` r
library(recalc)

res <- recalc_independent_t(
  m1 = 10.3, m2 = 8.7,
  sd1 = 3.1, sd2 = 2.8,
  n1 = 50,   n2 = 48,
  m_digits  = 1,
  sd_digits = 1,
  d = 0.55, d_ci_lower = 0.20, d_ci_upper = 0.90, d_digits = 2,
  p = 0.0024, p_digits = 4
)

res$reproduced          # p_inbounds / d_inbounds verdicts + interval bounds
plot_multiverse_p(res)  # every recalculated p against the reported one
plot_multiverse_d(res)
```

Setting `include_se_sd_confusion = TRUE` adds a branch to the multiverse in which the reported SDs are reinterpreted as SEs — a common source of irreproducible *t*-tests.

### Contingency tables

``` r
counts <- matrix(c(30, 20,
                   15, 35), nrow = 2, byrow = TRUE)

recalc_chisq(counts = counts, p = 0.002, p_digits = 3)$reproduced
```

## Function reference

| function | purpose |
|---|---|
| **Independent-samples *t*-test** | |
| `recalc_independent_t()` | combined multiverse of *p* and *d* for an independent-samples *t*-test |
| `recalc_independent_t_p()` | *p*-value multiverse only |
| `recalc_independent_t_d()` | Cohen's *d* / Hedges' *g* and CI multiverse only |
| `plot_multiverse_p()`, `plot_multiverse_d()` | plot the *p* / *d* multiverse against the reported value |
| **Contingency tables** | |
| `recalc_chisq()` / `recalc_chisq_p()` | multiverse of chi-squared / association-test *p*-values (Pearson, Yates, likelihood-ratio, Fisher, …) |
| **Regression** | |
| `recalc_regression_from_cor_mat()` | refit a linear regression from a reported correlation matrix |
| `audit_regression()`, `summarise_audit()` | run every applicable regression check and summarise it |
| `recalc_r2_from_f()`, `recalc_f_from_r2()`, `recalc_adj_r2()`, `recalc_delta_f()` | convert between R², *F*, adjusted R², and ΔF |
| `recalc_r2_from_betas_corrs()`, `recalc_r2_from_blocks()`, `recalc_r2_floor_from_correlations()` | recompute / bound R² from betas, correlations, or hierarchical blocks |
| `recalc_beta_from_b()`, `recalc_t_from_b_se()`, `recalc_ci_from_b_se()`, `recalc_se_from_ci()`, `recalc_p_from_t_df()` | coefficient-level conversions |
| `recalc_semipartial_r2_from_t()`, `recalc_partial_eta_from_f()`, `recalc_p_from_f()`, `recalc_p_model_from_r2()` | effect-size and model-level recalculations |
| `diagnose_beta_label()`, `diagnose_r2_label()`, `diagnose_p_tails()` | flag mislabeled coefficients / R² / tail conventions |
| **Mediation** | |
| `recalc_mediation()`, `recalc_mediation_ab()`, `recalc_mediation_p()` | recalculate an indirect effect and its Sobel-family *p*-value |
| **Reliability** | |
| `recalc_alpha_from_cor_sd()` | recalculate Cronbach's alpha from a correlation matrix and SDs |
| **Within-person / pre-post** | |
| `recalc_prepost_r()`, `recalc_prepost_r_from_f()`, `recalc_change_sd_from_r()` | implied pre-post correlation / change-score SD |
| **Subgroups** | |
| `recalc_total_from_subgroups()`, `recalc_missing_subgroup()` | recover overall or missing-subgroup N / M / SD |
| `recalc_sd_concentration()` | test for insufficient variance among subgroup SDs (exploratory) |
| **Helpers** | |
| `triangle_to_cor_matrix()` | expand a correlation triangle to a full symmetric matrix |

## Related projects

- Mark Bolland's [`reappraised`](https://reappraised.shinyapps.io/check_p_vals_cont/) R package and Shiny app aim to do the same as this package, with different degrees of analytic flexibility.
- Different versions of between-groups Cohen's *d*: <https://rpubs.com/metinbulus/welch>
- Different versions of within-subjects Cohen's *d*: <https://github.com/ianhussey/versions-of-cohens-d>
- Extracting pre-post correlations from summary statistics <https://matthewbjane.quarto.pub/pre-post-correlations/>
- Lisa DeBruine's [`faux`](https://github.com/debruine/faux) package has some error-detection functions.
- Lisa DeBruine's [`within`](https://github.com/debruine/within) repo assesses the plausibility of a within-subject *t*-test by showing the range of possible between-timepoint correlations given the reported results.

## Suggested citation

Hussey, I. (2026). *recalc*: Trustworthiness assessments for summary statistics using recalculation checks. [Computer software]. <https://github.com/ianhussey/recalc>
[doi:10.5281/zenodo.21439896](https://doi.org/10.5281/zenodo.21439896)

## License

Code is MIT licensed © Ian Hussey (2024).

Text and images are CC BY 4.0.
