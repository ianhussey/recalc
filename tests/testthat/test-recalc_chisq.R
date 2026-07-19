test_that("recalc_chisq_p reproduces base-R tests for a 2x2", {
  tab <- matrix(c(3, 12, 9, 6), nrow = 2, byrow = TRUE)
  r <- recalc_chisq_p(counts = tab)
  got <- function(m) r$p_results$p_unrounded[r$p_results$method == m]

  expect_equal(
    got("pearson"),
    suppressWarnings(chisq.test(tab, correct = FALSE)$p.value)
  )
  expect_equal(
    got("yates"),
    suppressWarnings(chisq.test(tab, correct = TRUE)$p.value)
  )
  expect_equal(got("fisher"), fisher.test(tab)$p.value)
})

test_that("n1_chisq is (N-1)/N * Pearson", {
  tab <- matrix(c(3, 12, 9, 6), nrow = 2, byrow = TRUE)
  N <- sum(tab)
  pear_stat <- suppressWarnings(chisq.test(tab, correct = FALSE)$statistic)
  expect_equal(
    recalc_chisq_p(counts = tab, p_methods = "n1_chisq")$p_results$p_unrounded[
      1
    ],
    unname(pchisq((N - 1) / N * pear_stat, 1, lower.tail = FALSE))
  )
})

test_that("fisher_central is the equal-tailed (2*min one-sided) convention", {
  tab <- matrix(c(1, 8, 15, 6), nrow = 2, byrow = TRUE) # asymmetric: central != minlike
  a <- tab[1, 1]
  r1 <- sum(tab[1, ])
  r2 <- sum(tab[2, ])
  c1 <- sum(tab[, 1])
  central <- min(
    1,
    2 *
      min(phyper(a, r1, r2, c1), phyper(a - 1, r1, r2, c1, lower.tail = FALSE))
  )
  expect_equal(
    recalc_chisq_p(
      counts = tab,
      p_methods = "fisher_central"
    )$p_results$p_unrounded[1],
    central
  )
  # on this asymmetric table it must differ from the minlike "fisher"
  expect_false(isTRUE(all.equal(
    recalc_chisq_p(counts = tab, p_methods = "fisher")$p_results$p_unrounded[1],
    central
  )))
})

test_that("wald_z is the unpooled two-proportion z (and != pooled/Pearson)", {
  tab <- matrix(c(3, 12, 9, 6), nrow = 2, byrow = TRUE)
  n1 <- sum(tab[1, ])
  n2 <- sum(tab[2, ])
  p1 <- tab[1, 1] / n1
  p2 <- tab[2, 1] / n2
  z <- (p1 - p2) / sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
  expect_equal(
    recalc_chisq_p(counts = tab, p_methods = "wald_z")$p_results$p_unrounded[1],
    2 * pnorm(-abs(z))
  )
  expect_false(isTRUE(all.equal(
    recalc_chisq_p(counts = tab, p_methods = "wald_z")$p_results$p_unrounded[1],
    suppressWarnings(chisq.test(tab, correct = FALSE)$p.value)
  )))
})

test_that("the default method set spans all nine methods, dropping 2x2-only for r x c", {
  r2 <- recalc_chisq_p(counts = matrix(c(3, 12, 9, 6), 2, byrow = TRUE))
  expect_setequal(
    unique(r2$p_results$method),
    c(
      "pearson",
      "yates",
      "likelihood_ratio",
      "n1_chisq",
      "fisher",
      "fisher_central",
      "fisher_midp",
      "fisher_midp_sas",
      "wald_z"
    )
  )
  r3 <- suppressMessages(
    recalc_chisq_p(counts = matrix(c(10, 20, 15, 12, 8, 18), 3, byrow = TRUE))
  )
  # 2x2-only methods gone; scaled-Pearson n1_chisq retained
  expect_true(all(
    c("pearson", "likelihood_ratio", "fisher", "n1_chisq") %in%
      r3$p_results$method
  ))
  expect_false(any(
    c(
      "yates",
      "fisher_central",
      "fisher_midp",
      "fisher_midp_sas",
      "wald_z"
    ) %in%
      r3$p_results$method
  ))
})

test_that("new methods compose into the gap-aware union reproduced summary", {
  r <- recalc_chisq_p(
    counts = matrix(c(3, 12, 9, 6), 2, byrow = TRUE),
    p = 0.05,
    p_digits = 3
  )
  expect_true(all(
    c(
      "p_inbounds",
      "p_inbounds_hull",
      "covered_width",
      "hull_width",
      "covered_fraction"
    ) %in%
      names(r$reproduced)
  ))
  # union is a subset of the hull, so never more permissive
  if (!is.na(r$reproduced$p_inbounds) && r$reproduced$p_inbounds) {
    expect_true(r$reproduced$p_inbounds_hull)
  }
})
