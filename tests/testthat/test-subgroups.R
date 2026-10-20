test_that("propagate_intervals matches brute-force corner enumeration", {
  inputs <- list(
    a = c(1, 2),
    b = c(3, 3), # degenerate interval: an exact value
    c = c(-1, 0.5)
  )
  fn <- function(a, b, c) a * b - c^2
  grid <- expand.grid(inputs)
  brute <- vapply(
    seq_len(nrow(grid)),
    function(i) fn(grid$a[i], grid$b[i], grid$c[i]),
    numeric(1)
  )
  expect_equal(
    propagate_intervals(fn, inputs),
    c(lower = min(brute), upper = max(brute))
  )
})

test_that("exact counts flag an N that rounding of the ns would hide", {
  args <- list(
    subgroup_ns = c(16, 16),
    subgroup_means = c(4.50, 4.90),
    overall_n = 33,
    mean_digits = 2
  )
  exact <- do.call(recalc_total_from_subgroups, c(args, n_digits = Inf))
  rounded <- do.call(recalc_total_from_subgroups, c(args, n_digits = 0))

  expect_equal(exact$recalculated_lower[1], 32)
  expect_equal(exact$recalculated_upper[1], 32)
  expect_false(exact$consistent[1])
  # n_digits = 0 gives each n a +/-0.5 interval, so 16 + 16 can reach 33
  expect_true(rounded$consistent[1])
})

test_that("propagate_intervals can leave out undefined corners", {
  inputs <- list(a = c(1, 2), b = c(3, 4))
  fn <- function(a, b) if (a == 1) NA_real_ else a * b
  expect_equal(
    propagate_intervals(fn, inputs, na.rm = TRUE),
    c(lower = 6, upper = 8)
  )
  expect_equal(
    propagate_intervals(function(a, b) NA_real_, inputs, na.rm = TRUE),
    c(lower = NA_real_, upper = NA_real_)
  )
})

missing_sd_row <- function(overall_sd, ...) {
  suppressWarnings(recalc_missing_subgroup(
    reported_ns = c(71, 62),
    reported_means = c(4.58, 4.03),
    reported_sds = c(1.45, 1.19),
    overall_n = 139,
    overall_mean = 4.30,
    overall_sd = overall_sd,
    n_digits = Inf,
    mean_digits = 2,
    sd_digits = 2,
    ...
  ))[3, ]
}

test_that("an SD feasible at only some roundings is not called impossible", {
  # at overall SD = 1.34, the implied variance is negative at some corners of
  # the rounding box but positive at most of them
  res <- missing_sd_row(1.34)
  expect_equal(res$recalculated_lower, 0)
  expect_equal(res$recalculated_upper, 1.2036, tolerance = 1e-4)
  expect_true(res$consistent)
  expect_true(missing_sd_row(1.34, scale_min = 1, scale_max = 7)$consistent)
})

test_that("an SD infeasible at every rounding is still impossible", {
  for (sd in c(1.31, 0.10)) {
    res <- missing_sd_row(sd)
    expect_true(is.nan(res$recalculated_lower))
    expect_true(is.nan(res$recalculated_upper))
    expect_false(res$consistent)
  }
  # and the Bhatia-Davis bound still applies
  expect_false(missing_sd_row(1.75, scale_min = 1, scale_max = 7)$consistent)
})

test_that("a missing group of one participant is checked without an SD", {
  a <- rep(1:7, length.out = 40)
  b <- rep(c(2, 3, 5, 6), length.out = 30)
  everyone <- c(a, b, 7)
  check <- function(sd_shift) {
    recalc_missing_subgroup(
      reported_ns = c(40, 30),
      reported_means = round(c(mean(a), mean(b)), 2),
      reported_sds = round(c(sd(a), sd(b)), 2),
      overall_n = 71,
      overall_mean = round(mean(everyone), 2),
      overall_sd = round(sd(everyone), 2) + sd_shift,
      n_digits = Inf,
      mean_digits = 2,
      sd_digits = 2
    )[3, ]
  }
  real <- check(0)
  expect_true(is.na(real$recalculated_lower))
  expect_true(real$consistent)
  expect_false(check(0.05)$consistent)
})

test_that("exact counts give an exact implied missing-subgroup n", {
  res <- recalc_missing_subgroup(
    reported_ns = 40,
    reported_means = 4.71,
    reported_sds = 2.04,
    overall_n = 80,
    overall_mean = 4.96,
    overall_sd = 2.18,
    n_digits = Inf,
    mean_digits = 2,
    sd_digits = 2
  )
  expect_equal(res$recalculated_lower[1], 40)
  expect_equal(res$recalculated_upper[1], 40)
  expect_true(all(res$consistent[c(1, 3)]))
})
