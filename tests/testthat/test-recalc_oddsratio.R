helper_tab <- function() {
  matrix(c(20, 30, 10, 40), nrow = 2, byrow = TRUE) # a=20,b=30,c=10,d=40
}

or_row <- function(res, estimator, ci_method, orientation = "as_is") {
  d <- res$or_results
  d[
    d$estimator == estimator &
      d$ci_method == ci_method &
      d$orientation == orientation,
    ,
    drop = FALSE
  ]
}

# --- point estimates ---------------------------------------------------------

test_that("unconditional MLE is ad/bc and conditional MLE matches fisher.test", {
  tab <- helper_tab()
  res <- recalc_oddsratio_p(counts = tab, or_orientation = "as_is")
  expect_equal(
    or_row(res, "uncond_mle", "wald_logor")$or_unrounded,
    (20 * 40) / (30 * 10)
  )
  expect_equal(
    or_row(res, "cond_mle", "wald_logor")$or_unrounded,
    unname(fisher.test(tab)$estimate)
  )
})

# --- CIs vs gold-standard implementations ------------------------------------

test_that("Cornfield CI equals fisher.test conf.int", {
  tab <- helper_tab()
  res <- recalc_oddsratio_p(
    counts = tab,
    or_estimators = "uncond_mle",
    or_orientation = "as_is"
  )
  r <- or_row(res, "uncond_mle", "cornfield_exact")
  ci <- fisher.test(tab)$conf.int
  expect_equal(r$ci_lower_unrounded, ci[1])
  expect_equal(r$ci_upper_unrounded, ci[2])
})

test_that("Wald CI and p equal the closed-form log-OR normal approximation", {
  tab <- helper_tab()
  logor <- log((20 * 40) / (30 * 10))
  se <- sqrt(1 / 20 + 1 / 30 + 1 / 10 + 1 / 40)
  res <- recalc_oddsratio_p(
    counts = tab,
    or_estimators = "uncond_mle",
    or_orientation = "as_is"
  )
  r <- or_row(res, "uncond_mle", "wald_logor")
  expect_equal(r$ci_lower_unrounded, exp(logor - qnorm(0.975) * se))
  expect_equal(r$ci_upper_unrounded, exp(logor + qnorm(0.975) * se))
  expect_equal(r$p_unrounded, 2 * pnorm(-abs(logor / se)))
})

test_that("profile-likelihood CI equals confint() on a logistic glm", {
  tab <- helper_tab()
  fit <- glm(
    cbind(c(20, 10), c(30, 40)) ~ factor(c("g1", "g2")),
    family = binomial()
  )
  ci_coef <- suppressMessages(confint(fit, 2L))
  expected <- exp(-rev(ci_coef)) # grp coef is -log(OR_native)
  res <- recalc_oddsratio_p(
    counts = tab,
    or_estimators = "uncond_mle",
    or_orientation = "as_is"
  )
  r <- or_row(res, "uncond_mle", "profile_likelihood")
  expect_equal(r$ci_lower_unrounded, unname(expected[1]), tolerance = 1e-4)
  expect_equal(r$ci_upper_unrounded, unname(expected[2]), tolerance = 1e-4)
})

test_that("mid-p CI sits inside the Cornfield exact CI and contains the OR", {
  tab <- helper_tab()
  res <- recalc_oddsratio_p(
    counts = tab,
    or_estimators = "uncond_mle",
    or_orientation = "as_is"
  )
  midp <- or_row(res, "uncond_mle", "midp_exact")
  corn <- or_row(res, "uncond_mle", "cornfield_exact")
  or <- midp$or_unrounded
  expect_gt(midp$ci_lower_unrounded, corn$ci_lower_unrounded)
  expect_lt(midp$ci_upper_unrounded, corn$ci_upper_unrounded)
  expect_lt(midp$ci_lower_unrounded, or)
  expect_gt(midp$ci_upper_unrounded, or)
})

# --- orientation -------------------------------------------------------------

test_that("orientation 'both' carries the reciprocal; 'as_is' does not", {
  tab <- helper_tab()
  both <- recalc_oddsratio_p(counts = tab, or_orientation = "both")
  asis <- recalc_oddsratio_p(counts = tab, or_orientation = "as_is")
  expect_true(all(c("as_is", "reciprocal") %in% both$or_results$orientation))
  expect_setequal(asis$or_results$orientation, "as_is")

  recip <- or_row(both, "uncond_mle", "wald_logor", "reciprocal")
  native <- or_row(both, "uncond_mle", "wald_logor", "as_is")
  expect_equal(recip$or_unrounded, 1 / native$or_unrounded)
  # reciprocal CI bounds are the inverted-and-swapped native bounds
  expect_equal(recip$ci_lower_unrounded, 1 / native$ci_upper_unrounded)
  expect_equal(recip$ci_upper_unrounded, 1 / native$ci_lower_unrounded)
  # Wald p is orientation-symmetric
  expect_equal(recip$p_unrounded, native$p_unrounded)
})

# --- continuity correction ---------------------------------------------------

test_that("continuity correction follows the on_zero/always/never policy", {
  z <- matrix(c(0, 10, 8, 12), nrow = 2, byrow = TRUE)
  on_zero <- recalc_chisq_p(
    counts = z,
    include_or = TRUE,
    or_estimators = "uncond_mle",
    or_ci_methods = "wald_logor",
    or_orientation = "as_is"
  )
  expect_true(all(on_zero$or_results$cc_label == "cc:+0.5"))
  expect_true(is.finite(on_zero$or_results$or_unrounded[1]))

  never <- recalc_chisq_p(
    counts = z,
    include_or = TRUE,
    or_estimators = "uncond_mle",
    or_ci_methods = "wald_logor",
    or_orientation = "as_is",
    or_cc_apply = "never"
  )
  expect_true(all(never$or_results$cc_label == "cc:none"))
  expect_true(is.na(never$or_results$p_unrounded[1]))

  # No zero cell + on_zero => no correction applied.
  nz <- recalc_chisq_p(
    counts = helper_tab(),
    include_or = TRUE,
    or_estimators = "uncond_mle",
    or_ci_methods = "wald_logor",
    or_orientation = "as_is"
  )
  expect_true(all(nz$or_results$cc_label == "cc:none"))

  # always with cc = 1 adds 1 to every cell.
  one <- recalc_chisq_p(
    counts = helper_tab(),
    include_or = TRUE,
    or_estimators = "uncond_mle",
    or_ci_methods = "wald_logor",
    or_orientation = "as_is",
    or_cc = 1,
    or_cc_apply = "always"
  )
  expect_true(all(one$or_results$cc_label == "cc:+1"))
  expect_equal(one$or_results$or_unrounded[1], (21 * 41) / (31 * 11))
})

# --- opt-in / invariance -----------------------------------------------------

test_that("the OR track is off by default and leaves the p output untouched", {
  tab <- helper_tab()
  off <- recalc_chisq_p(counts = tab, p = 0.02, p_digits = 2)
  expect_setequal(names(off), c("reproduced", "p_results", "method_intervals"))
  expect_false(any(grepl("^or_|odds_ratio", names(off$reproduced))))

  on <- recalc_chisq_p(counts = tab, p = 0.02, p_digits = 2, include_or = TRUE)
  expect_equal(on$p_results, off$p_results)
  expect_true(all(c("or_results", "or_method_intervals") %in% names(on)))
})

test_that("supplying a reported OR auto-activates the track", {
  res <- recalc_chisq_p(counts = helper_tab(), odds_ratio = 2.67, or_digits = 2)
  expect_true("or_results" %in% names(res))
  expect_equal(res$reproduced$odds_ratio, 2.67)
})

test_that("statistic and df are surfaced in method_intervals", {
  res <- recalc_chisq_p(counts = helper_tab())
  mi <- res$method_intervals
  expect_true(all(c("statistic_min", "statistic_max", "df") %in% names(mi)))
  expect_equal(mi$df[mi$method == "pearson"], 1)
  expect_true(is.na(mi$statistic_min[mi$method == "fisher"]))
})

# --- reproduced decisions ----------------------------------------------------

test_that("union OR decision covers both native and reciprocal clusters", {
  tab <- helper_tab()
  native <- (20 * 40) / (30 * 10)
  res <- recalc_oddsratio_p(counts = tab, odds_ratio = round(native, 2))
  expect_true(res$reproduced$or_inbounds)
  # reported as 1/OR is still reproducible under orientation 'both'
  res2 <- recalc_oddsratio_p(counts = tab, odds_ratio = round(1 / native, 2))
  expect_true(res2$reproduced$or_inbounds)
})

test_that("or_p_inbounds reflects the Wald-log-OR p, not the chi-square p", {
  tab <- helper_tab()
  logor <- log((20 * 40) / (30 * 10))
  se <- sqrt(1 / 20 + 1 / 30 + 1 / 10 + 1 / 40)
  wald_p <- 2 * pnorm(-abs(logor / se))
  res <- recalc_oddsratio_p(
    counts = tab,
    p = round(wald_p, 3),
    p_digits = 3
  )
  expect_true(res$reproduced$or_p_inbounds)
  # a p the Wald test cannot produce is out of the OR-p bounds
  res2 <- recalc_oddsratio_p(counts = tab, p = 0.001, p_digits = 3)
  expect_false(res2$reproduced$or_p_inbounds)
})

# --- guards / wrapper --------------------------------------------------------

test_that("OR track is skipped with a message for tables larger than 2x2", {
  big <- matrix(c(10, 20, 30, 15, 25, 5), nrow = 2, byrow = TRUE)
  expect_message(
    res <- recalc_chisq_p(counts = big, include_or = TRUE),
    "only defined for a 2x2"
  )
  expect_false("or_results" %in% names(res))
})

test_that("recalc_oddsratio_p equals include_or = TRUE and orders OR first", {
  tab <- helper_tab()
  w <- recalc_oddsratio_p(counts = tab)
  e <- recalc_chisq_p(counts = tab, include_or = TRUE)
  expect_equal(w$or_results, e$or_results)
  expect_equal(names(w$reproduced)[1], "odds_ratio")
})

test_that("invalid OR method names error", {
  expect_error(
    recalc_chisq_p(counts = helper_tab(), include_or = TRUE, or_estimators = "bogus"),
    "or_estimators"
  )
  expect_error(
    recalc_chisq_p(counts = helper_tab(), include_or = TRUE, or_ci_methods = "bogus"),
    "or_ci_methods"
  )
})

# --- proportion route --------------------------------------------------------

test_that("the OR track works through the proportion-reconstruction route", {
  res <- recalc_oddsratio_p(
    prop1 = 0.40,
    prop2 = 0.20,
    n1 = 50,
    n2 = 50,
    prop_digits = 2,
    or_orientation = "as_is"
  )
  expect_true("or_results" %in% names(res))
  expect_true(nrow(res$or_results) > 0)
  expect_true(all(is.finite(res$or_results$or_unrounded)))
})

test_that("profile-likelihood CI has a finite limit and an infinite one with a zero cell", {
  m <- function(v) matrix(v, 2, byrow = TRUE)
  ci <- .or_ci_profile(m(c(4, 11, 0, 15)), 0.05)
  # the OR is infinite, so the upper limit is infinite (NA), and the lower
  # limit is where the likelihood-ratio deviance reaches the chi-square quantile
  expect_true(is.na(ci[2]))
  expect_equal(ci[1], 1.972033, tolerance = 1e-5)
  b0 <- seq(-40, 10, by = 1e-4)
  ll <- 4 * (b0 + log(ci[1])) - 15 * log1p(exp(b0 + log(ci[1]))) -
    15 * log1p(exp(b0))
  ll_max <- 4 * log(4 / 15) + 11 * log(11 / 15)
  expect_equal(2 * (ll_max - max(ll)), qchisq(0.95, 1), tolerance = 1e-6)

  # with the rows swapped, the CI is the reciprocal
  swapped <- .or_ci_profile(m(c(0, 15, 4, 11)), 0.05)
  expect_equal(swapped[1], 0)
  expect_equal(swapped[2], 1 / ci[1], tolerance = 1e-6)

  # an empty column has no OR at all
  expect_true(all(is.na(.or_ci_profile(m(c(0, 5, 0, 5)), 0.05))))

  # through the multiverse
  res <- recalc_oddsratio(counts = c(4, 11, 0, 15), n_rows = 2, n_cols = 2)
  prof <- res$or_results[
    res$or_results$ci_method == "profile_likelihood" &
      res$or_results$orientation == "as_is",
  ]
  expect_true(all(abs(prof$ci_lower_unrounded - 1.972033) < 1e-5))
})

test_that("direct LR inversion agrees with confint() on tables without zeros", {
  m <- function(v) matrix(v, 2, byrow = TRUE)
  for (v in list(c(13, 123, 2, 135), c(10, 20, 15, 25), c(3, 7, 8, 2))) {
    expect_equal(
      .or_ci_profile_lr(m(v), 0.05),
      .or_ci_profile(m(v), 0.05),
      tolerance = 1e-3
    )
  }
})
