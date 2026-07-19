# recalc roadmap / TODO

Development notes moved out of the README.

- Expand to include the dependent (paired) t-test. The `within` repo has some
  code for this, but the basis should be the code in this repo. Add the choice
  of also using the reported t value, or p value, or pre-post correlation, or a
  reasonable estimate of it (e.g. .5 to .7, or perhaps .25 to .9) to reproduce
  the p and d. Produce estimates of both d_z and d_av.
- Use the independent t-test code to expand to one-way and two-way ANOVA,
  following Nick Brown's original goal.
- Document the regression recalc.
- Add regression-from-correlation-matrix recalc (done: see
  `recalc_regression_from_cor_mat()`).
- Add mediation analysis recalc, e.g. Rutherford et al. 2017 (done: see
  `recalc_mediation()`).
- **Integrate granularity and bounds tests more generally.** GRIM / GRIMMER-style
  granularity checks (does a reported mean / SD have a denominator that a stated
  N can produce?) and physical-bounds checks (is the reported mean within the
  scale's logical range; is the reported SD within the Bhatia-Davis bound implied
  by the scale endpoints and the reported mean / N) are currently only wired into
  `recalc_missing_subgroup()` (E2). They apply equally to E1's recalculated
  quantities, to `recalc_independent_t()`'s reported group statistics, and to any
  other check whose inputs include means and SDs on a bounded scale. Factor
  `scale_min` / `scale_max` arguments, a shared Bhatia-Davis bound helper, and a
  granularity helper into a common layer used by every reported-descriptive check.
- Validate `recalc_sd_concentration()` (E3): the lower-tail Bartlett
  interpretation as a TIVA-style forensic signal is currently flagged as
  exploratory in its docstring and vignette. Run null-calibration sims at
  realistic (k, n_g) and under mild non-normality; characterize power against
  plausible fabricated-SD patterns; only then graduate it from "diagnostic" to
  "check" and consider inclusion in any `audit_*()` wrapper.
