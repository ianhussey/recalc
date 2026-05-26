# Tests for recalc_independent_t_p, recalc_independent_t_d, and the
# combined wrapper.
#
# Strategy: pick a fixed (m1, m2, sd1, sd2, n1, n2) reported at full
# precision (so the rounding grid collapses to the reported corner) and
# check that the multiverse's "reported|reported|reported|reported"
# row matches reference values from base R t.test() and known
# closed-form d / Hedges' g calculations. This isolates the analytic
# choices being exposed from the rounding multiverse.

# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

# A fixed scenario used across many tests. Choose values that have
# trailing zero so m_digits / sd_digits match the apparent precision.
.scn <- function() {
  list(m1 = 10.30, m2 = 8.71, sd1 = 3.12, sd2 = 2.80, n1 = 50, n2 = 48,
       m_digits = 2, sd_digits = 2)
}

# Pull the "reported corner" row(s) from a multiverse result frame.
.reported_rows <- function(df) {
  df[df$input_adj_stats == "m1:reported|m2:reported|sd1:reported|sd2:reported", ,
     drop = FALSE]
}

# ----------------------------------------------------------------------
# P-value reproduction against base R t.test
# ----------------------------------------------------------------------

test_that("Welch and Student two-sided p match base R t.test", {
  s <- .scn()
  res <- recalc_independent_t_p(
    m1 = s$m1, m2 = s$m2, sd1 = s$sd1, sd2 = s$sd2, n1 = s$n1, n2 = s$n2,
    m_digits = s$m_digits, sd_digits = s$sd_digits,
    p_methods = c("student_t", "welch_t"),
    alternative = "two.sided"
  )

  rep <- .reported_rows(res$p_results)
  p_student <- rep$p_unrounded[rep$p_method == "student_t"]
  p_welch   <- rep$p_unrounded[rep$p_method == "welch_t"]

  # Simulate equivalent data to drive base R's t.test from summary stats:
  # easier — replicate the t-statistic + df closed form.
  sp <- sqrt(((s$n1 - 1) * s$sd1^2 + (s$n2 - 1) * s$sd2^2) /
               (s$n1 + s$n2 - 2))
  t_s <- (s$m1 - s$m2) / (sp * sqrt(1 / s$n1 + 1 / s$n2))
  p_s_ref <- 2 * pt(-abs(t_s), df = s$n1 + s$n2 - 2)

  v1 <- s$sd1^2; v2 <- s$sd2^2
  se_w <- sqrt(v1 / s$n1 + v2 / s$n2)
  t_w <- (s$m1 - s$m2) / se_w
  df_w <- (v1 / s$n1 + v2 / s$n2)^2 /
    ((v1^2 / (s$n1^2 * (s$n1 - 1))) + (v2^2 / (s$n2^2 * (s$n2 - 1))))
  p_w_ref <- 2 * pt(-abs(t_w), df = df_w)

  expect_equal(p_student, p_s_ref, tolerance = 1e-12)
  expect_equal(p_welch,   p_w_ref, tolerance = 1e-12)
})

test_that("one-sided alternatives match the directed tail", {
  s <- .scn()
  res <- recalc_independent_t_p(
    m1 = s$m1, m2 = s$m2, sd1 = s$sd1, sd2 = s$sd2, n1 = s$n1, n2 = s$n2,
    m_digits = s$m_digits, sd_digits = s$sd_digits,
    p_methods = "student_t",
    alternative = c("two.sided", "less", "greater")
  )
  rep <- .reported_rows(res$p_results)

  sp <- sqrt(((s$n1 - 1) * s$sd1^2 + (s$n2 - 1) * s$sd2^2) /
               (s$n1 + s$n2 - 2))
  t_s <- (s$m1 - s$m2) / (sp * sqrt(1 / s$n1 + 1 / s$n2))
  df  <- s$n1 + s$n2 - 2

  expect_equal(rep$p_unrounded[rep$alternative == "two.sided"],
               2 * pt(-abs(t_s), df = df), tolerance = 1e-12)
  expect_equal(rep$p_unrounded[rep$alternative == "less"],
               pt(t_s, df = df, lower.tail = TRUE), tolerance = 1e-12)
  expect_equal(rep$p_unrounded[rep$alternative == "greater"],
               pt(t_s, df = df, lower.tail = FALSE), tolerance = 1e-12)

  # The one-sided p-values should sum to ~1
  expect_equal(rep$p_unrounded[rep$alternative == "less"] +
                 rep$p_unrounded[rep$alternative == "greater"],
               1, tolerance = 1e-12)
})

test_that("tail precision: very small two-sided p does not underflow to 0", {
  # Large t -> tiny p. The old form 2*(1 - pt(abs(t), df)) underflows;
  # the new form 2*pt(-abs(t), df) keeps precision.
  s <- list(m1 = 10, m2 = 1, sd1 = 1, sd2 = 1, n1 = 50, n2 = 50,
            m_digits = 0, sd_digits = 0)
  res <- recalc_independent_t_p(
    m1 = s$m1, m2 = s$m2, sd1 = s$sd1, sd2 = s$sd2, n1 = s$n1, n2 = s$n2,
    m_digits = s$m_digits, sd_digits = s$sd_digits,
    p_methods = "student_t", alternative = "two.sided"
  )
  rep <- .reported_rows(res$p_results)
  p_val <- rep$p_unrounded[1]
  expect_true(p_val > 0)
  expect_true(p_val < 1e-30)
})

# ----------------------------------------------------------------------
# d / SE / CI variants
# ----------------------------------------------------------------------

test_that("d_formulas produce the documented denominators", {
  s <- .scn()
  res <- recalc_independent_t_d(
    m1 = s$m1, m2 = s$m2, sd1 = s$sd1, sd2 = s$sd2, n1 = s$n1, n2 = s$n2,
    m_digits = s$m_digits, sd_digits = s$sd_digits,
    d_formulas = c("pooled_df", "pooled_n", "simple_avg"),
    hedges_correction = "none",
    ci_methods = "wald_t",
    se_formulas = "cumming"
  )
  rep <- .reported_rows(res$d_results)

  sp_df  <- sqrt(((s$n1 - 1) * s$sd1^2 + (s$n2 - 1) * s$sd2^2) /
                   (s$n1 + s$n2 - 2))
  sp_n   <- sqrt((s$n1 * s$sd1^2 + s$n2 * s$sd2^2) / (s$n1 + s$n2))
  sp_avg <- sqrt((s$sd1^2 + s$sd2^2) / 2)

  expect_equal(unique(rep$d_unrounded[rep$d_formula == "pooled_df"]),
               (s$m1 - s$m2) / sp_df,  tolerance = 1e-12)
  expect_equal(unique(rep$d_unrounded[rep$d_formula == "pooled_n"]),
               (s$m1 - s$m2) / sp_n,   tolerance = 1e-12)
  expect_equal(unique(rep$d_unrounded[rep$d_formula == "simple_avg"]),
               (s$m1 - s$m2) / sp_avg, tolerance = 1e-12)
})

test_that("Hedges' correction variants match expected J values", {
  s <- .scn()
  res <- recalc_independent_t_d(
    m1 = s$m1, m2 = s$m2, sd1 = s$sd1, sd2 = s$sd2, n1 = s$n1, n2 = s$n2,
    m_digits = s$m_digits, sd_digits = s$sd_digits,
    d_formulas = "pooled_df",
    hedges_correction = c("none", "approx", "exact"),
    ci_methods = "wald_t", se_formulas = "cumming"
  )
  rep <- .reported_rows(res$d_results)

  sp <- sqrt(((s$n1 - 1) * s$sd1^2 + (s$n2 - 1) * s$sd2^2) /
               (s$n1 + s$n2 - 2))
  d_raw <- (s$m1 - s$m2) / sp
  df <- s$n1 + s$n2 - 2
  J_approx <- 1 - 3 / (4 * df - 1)
  J_exact  <- exp(lgamma(df / 2) - 0.5 * log(df / 2) - lgamma((df - 1) / 2))

  expect_equal(unique(rep$d_unrounded[rep$hedges_correction == "none"]),
               d_raw, tolerance = 1e-12)
  expect_equal(unique(rep$d_unrounded[rep$hedges_correction == "approx"]),
               J_approx * d_raw, tolerance = 1e-12)
  expect_equal(unique(rep$d_unrounded[rep$hedges_correction == "exact"]),
               J_exact * d_raw, tolerance = 1e-12)

  # Approx and exact should be very close at df = 96.
  expect_lt(abs(J_approx - J_exact), 1e-4)
})

test_that("SE-of-d formulas differ for wald_t and match closed forms", {
  s <- .scn()
  res <- recalc_independent_t_d(
    m1 = s$m1, m2 = s$m2, sd1 = s$sd1, sd2 = s$sd2, n1 = s$n1, n2 = s$n2,
    m_digits = s$m_digits, sd_digits = s$sd_digits,
    d_formulas = "pooled_df", hedges_correction = "none",
    ci_methods = "wald_t",
    se_formulas = c("hedges_olkin", "cumming", "uncorrected")
  )
  rep <- .reported_rows(res$d_results)

  N <- s$n1 + s$n2
  sp <- sqrt(((s$n1 - 1) * s$sd1^2 + (s$n2 - 1) * s$sd2^2) /
               (s$n1 + s$n2 - 2))
  d  <- (s$m1 - s$m2) / sp
  df <- N - 2
  crit_t <- qt(0.975, df = df)

  se_ho <- sqrt(N / (s$n1 * s$n2) + d^2 / (2 * N))
  se_cu <- sqrt(N / (s$n1 * s$n2) + d^2 / (2 * (N - 2)))
  se_un <- sqrt(N / (s$n1 * s$n2))

  expect_equal(
    rep$ci_lower_unrounded[rep$se_formula == "hedges_olkin"][1],
    d - crit_t * se_ho, tolerance = 1e-12
  )
  expect_equal(
    rep$ci_lower_unrounded[rep$se_formula == "cumming"][1],
    d - crit_t * se_cu, tolerance = 1e-12
  )
  expect_equal(
    rep$ci_lower_unrounded[rep$se_formula == "uncorrected"][1],
    d - crit_t * se_un, tolerance = 1e-12
  )
})

test_that("wald_z uses z critical (not t)", {
  s <- .scn()
  res <- recalc_independent_t_d(
    m1 = s$m1, m2 = s$m2, sd1 = s$sd1, sd2 = s$sd2, n1 = s$n1, n2 = s$n2,
    m_digits = s$m_digits, sd_digits = s$sd_digits,
    d_formulas = "pooled_df", hedges_correction = "none",
    ci_methods = c("wald_t", "wald_z"),
    se_formulas = "cumming"
  )
  rep <- .reported_rows(res$d_results)

  N <- s$n1 + s$n2
  sp <- sqrt(((s$n1 - 1) * s$sd1^2 + (s$n2 - 1) * s$sd2^2) /
               (s$n1 + s$n2 - 2))
  d  <- (s$m1 - s$m2) / sp
  df <- N - 2
  se <- sqrt(N / (s$n1 * s$n2) + d^2 / (2 * (N - 2)))

  expect_equal(
    rep$ci_lower_unrounded[rep$ci_method == "wald_t"][1],
    d - qt(0.975, df = df) * se, tolerance = 1e-12
  )
  expect_equal(
    rep$ci_lower_unrounded[rep$ci_method == "wald_z"][1],
    d - qnorm(0.975) * se, tolerance = 1e-12
  )
})

test_that("NCT CI produces sensible bounds that bracket d", {
  s <- .scn()
  res <- recalc_independent_t_d(
    m1 = s$m1, m2 = s$m2, sd1 = s$sd1, sd2 = s$sd2, n1 = s$n1, n2 = s$n2,
    m_digits = s$m_digits, sd_digits = s$sd_digits,
    d_formulas = "pooled_df", hedges_correction = "none",
    ci_methods = "nct"
  )
  rep <- .reported_rows(res$d_results)
  d  <- rep$d_unrounded[1]
  lo <- rep$ci_lower_unrounded[1]
  hi <- rep$ci_upper_unrounded[1]
  expect_true(is.finite(lo) && is.finite(hi))
  expect_true(lo < d && d < hi)
})

# ----------------------------------------------------------------------
# Combined wrapper
# ----------------------------------------------------------------------

test_that("recalc_independent_t merges p and d sub-results", {
  s <- .scn()
  res <- recalc_independent_t(
    m1 = s$m1, m2 = s$m2, sd1 = s$sd1, sd2 = s$sd2, n1 = s$n1, n2 = s$n2,
    m_digits = s$m_digits, sd_digits = s$sd_digits,
    p = 0.01, p_digits = 2, p_operator = "less_than",
    d = 0.53, d_digits = 2
  )
  expect_true(all(c("reproduced", "d_results", "p_results") %in% names(res)))
  expect_true(nrow(res$reproduced) == 1L)
  expect_true(all(c("d", "min_d", "max_d", "d_inbounds",
                    "p_operator", "p", "min_p", "max_p", "p_inbounds")
                  %in% names(res$reproduced)))
})

# ----------------------------------------------------------------------
# Validation errors
# ----------------------------------------------------------------------

test_that("missing m_digits / sd_digits errors", {
  s <- .scn()
  expect_error(
    recalc_independent_t_p(m1 = s$m1, m2 = s$m2, sd1 = s$sd1, sd2 = s$sd2,
                           n1 = s$n1, n2 = s$n2, sd_digits = 2),
    "m_digits"
  )
  expect_error(
    recalc_independent_t_d(m1 = s$m1, m2 = s$m2, sd1 = s$sd1, sd2 = s$sd2,
                           n1 = s$n1, n2 = s$n2, m_digits = 2),
    "sd_digits"
  )
})

test_that("invalid method / alternative names error", {
  s <- .scn()
  expect_error(
    recalc_independent_t_p(m1 = s$m1, m2 = s$m2, sd1 = s$sd1, sd2 = s$sd2,
                           n1 = s$n1, n2 = s$n2,
                           m_digits = s$m_digits, sd_digits = s$sd_digits,
                           p_methods = "nope"),
    "p_methods"
  )
  expect_error(
    recalc_independent_t_p(m1 = s$m1, m2 = s$m2, sd1 = s$sd1, sd2 = s$sd2,
                           n1 = s$n1, n2 = s$n2,
                           m_digits = s$m_digits, sd_digits = s$sd_digits,
                           alternative = "bad"),
    "alternative"
  )
})

test_that("hedges_correction accepts logical (back-compat)", {
  s <- .scn()
  res_T <- recalc_independent_t_d(
    m1 = s$m1, m2 = s$m2, sd1 = s$sd1, sd2 = s$sd2, n1 = s$n1, n2 = s$n2,
    m_digits = s$m_digits, sd_digits = s$sd_digits,
    hedges_correction = TRUE
  )
  expect_true(all(res_T$d_results$hedges_correction %in% c("approx", "exact")))

  res_F <- recalc_independent_t_d(
    m1 = s$m1, m2 = s$m2, sd1 = s$sd1, sd2 = s$sd2, n1 = s$n1, n2 = s$n2,
    m_digits = s$m_digits, sd_digits = s$sd_digits,
    hedges_correction = FALSE
  )
  expect_true(all(res_F$d_results$hedges_correction == "none"))
})
