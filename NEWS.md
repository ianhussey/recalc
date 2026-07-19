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
