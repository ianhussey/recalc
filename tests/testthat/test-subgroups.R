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
