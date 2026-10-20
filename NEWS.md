# recalc 0.6.4

## New features

* `recalc_independent_t_d()` and `recalc_independent_t()` now check a reported
  CI for d. They accepted `d_ci_lower` and `d_ci_upper` before but did not use
  them. A CI depends on every analytic choice, including the CI method and the
  SE formula, so for each combination of choices the range of each bound
  across input rounding is rounded to `d_digits`. `ci_lower_inbounds` and
  `ci_upper_inbounds` say whether each bound lies in that range for at least
  one combination, and `ci_inbounds` whether a single combination contains
  both. The per-combination ranges and verdicts are returned as
  `ci_method_intervals`. The t-test Shiny app shows the CI verdict.

## Improvements

* The d multiverse behind `recalc_independent_t_d()` and
  `recalc_independent_t()` is built once instead of one data frame per
  combination, making it about 8 times faster (e.g. 7,290 combinations:
  ~2.5 s to ~0.3 s). The results are identical.

# recalc 0.6.3

## Bug fixes

* `recalc_missing_subgroup()` no longer calls the implied missing-subgroup SD
  impossible when the implied variance is negative at only some corners of the
  rounding box. It propagated the SD, so a single negative-variance corner made
  both bounds `NaN` and the row `FALSE`, even when most roundings gave a valid
  SD (e.g. 52 of 64 corners in the test case). The variance is now propagated:
  the row is `FALSE`, with `NaN` bounds, only if the variance is negative at
  every corner, and an interval that straddles 0 gives an SD interval starting
  at 0.
* A missing group of exactly one participant, which has no SD, is now checked
  by whether the sum of squares left for it can be exactly 0, rather than by
  dividing by zero. Its SD bounds are `NA`.
* `propagate_intervals()` gains `na.rm`, to leave out corners where a quantity
  is undefined.

# recalc 0.6.2

## Improvements

* `propagate_intervals()`, which every rounding-aware check uses, drops
  duplicate corners (an input with a zero-width interval contributes one
  corner, not two) and no longer indexes the corner grid row by row. The
  results are unchanged; checks with many inputs are much faster (e.g. six
  subgroups with SDs in `recalc_total_from_subgroups()`: ~17 s to ~0.2 s).

## Documentation

* The subgroup checks (`recalc_total_from_subgroups()`,
  `recalc_missing_subgroup()`, `recalc_sd_concentration()`) now recommend
  `n_digits = Inf`, treating sample sizes as the exact counts they are.
  `n_digits = 0` gives each n a +/-0.5 interval, which lets a reported N that
  differs from the sum of the subgroup ns by up to (k + 1)/2 pass as
  consistent (e.g. 16 + 16 against a reported N of 33). The examples, the
  subgroup vignette, and the whole-vs-subgroup Shiny app (which gains an
  "N values are exact counts" option, on by default) use exact counts.

# recalc 0.6.1

## New features

* Odds-ratio multiverse for 2x2 tables, as an opt-in companion to the
  contingency-table p-value multiverse. `recalc_chisq_p()` gains an
  `include_or` switch (and auto-enables when a reported `odds_ratio` /
  `or_ci_lower` / `or_ci_upper` is supplied); `recalc_oddsratio_p()` /
  `recalc_oddsratio()` are thin wrappers that enable it and list the OR columns
  first. The track reconciles the odds-ratio point estimate (unconditional and
  conditional MLE), its CI (Wald-on-log-OR, Cornfield/exact-conditional,
  mid-p exact, and profile-likelihood), and the Wald-log-OR p-value, with
  continuity correction (`or_cc` / `or_cc_apply`, including the +1 variant) and
  an OR/`1/OR` orientation axis. The chi-square association p decision is left
  untouched, and the same gap-aware per-method union logic is used for the OR
  estimate and p.
* `recalc_chisq_p()` now surfaces the chi-square `statistic_min`,
  `statistic_max`, and `df` alongside each method's p interval in
  `method_intervals`.

# recalc 0.6.0

First public release, prepared for CRAN submission.

## Features

* Independent-samples *t*-test recalculation: `recalc_independent_t()`,
  `recalc_independent_t_p()`, and `recalc_independent_t_d()` return the
  multiverse of *p*-values and Cohen's *d* / Hedges' *g* (with CIs) compatible
  with reported means, SDs, and *N*s, taking input/output rounding, Student's
  vs Welch's *t*, one- vs two-sided tests, and SD/SE confusion into account.
* Contingency-table recalculation: `recalc_chisq()` / `recalc_chisq_p()` across
  Pearson, Yates, likelihood-ratio, Fisher, and related test variants.
* Regression checks: `recalc_regression_from_cor_mat()`, `audit_regression()`,
  `summarise_audit()`, and a family of coefficient-, effect-size-, and
  model-level recalculations (`recalc_r2_from_f()`, `recalc_beta_from_b()`,
  `recalc_semipartial_r2_from_t()`, and others), plus label-confusion
  diagnostics (`diagnose_beta_label()`, `diagnose_r2_label()`,
  `diagnose_p_tails()`).
* Mediation: `recalc_mediation()`, `recalc_mediation_ab()`,
  `recalc_mediation_p()`.
* Reliability: `recalc_alpha_from_cor_sd()`.
* Within-person / pre-post: `recalc_prepost_r()`, `recalc_prepost_r_from_f()`,
  `recalc_change_sd_from_r()`.
* Subgroup consistency: `recalc_total_from_subgroups()`,
  `recalc_missing_subgroup()`, and the exploratory `recalc_sd_concentration()`.
* Visualisation: `plot_multiverse_p()` and `plot_multiverse_d()`.

## Documentation

* Added the `vignette("recalc")` getting-started guide.
