test_that("recalc_prepost_r reproduces the identity (point, mid-interval)", {
  # r derived from SDs built to be mutually exact at a known r
  for (r in c(-0.6, -0.2, 0, 0.3, 0.7, 0.9)) {
    sd_pre <- 5; sd_post <- 8
    sd_chg <- sqrt(sd_pre^2 + sd_post^2 - 2 * r * sd_pre * sd_post)
    res <- recalc_prepost_r(sd_pre, sd_post, round(sd_chg, 6),
                            sd_pre_digits = 6, sd_post_digits = 6,
                            sd_change_digits = 6)
    mid <- (res$recalculated_lower + res$recalculated_upper) / 2
    expect_equal(mid, r, tolerance = 1e-4)
  }
})

test_that("recalc_prepost_r reproduces the Gauhar (2016) TSC-40 cases", {
  tab <- recalc_prepost_r(6.77, 7.45, 13.01,
                          sd_pre_digits = 2, sd_post_digits = 2, sd_change_digits = 2)
  mid_tab <- (tab$recalculated_lower + tab$recalculated_upper) / 2
  expect_equal(mid_tab, -0.673, tolerance = 0.005)   # table SDs -> negative r
  expect_true(tab$recalculated_lower < 0)

  txt <- recalc_prepost_r(14.58, 7.79, 13.01,
                          sd_pre_digits = 2, sd_post_digits = 2, sd_change_digits = 2)
  mid_txt <- (txt$recalculated_lower + txt$recalculated_upper) / 2
  expect_equal(mid_txt, 0.458, tolerance = 0.005)    # text SDs -> normal r
})

test_that("recalc_prepost_r flags impossibility via the interval (|r| > 1)", {
  # change SD below |sd_pre - sd_post| -> r > 1 throughout
  hi <- recalc_prepost_r(8, 5, 2.0,
                         sd_pre_digits = 1, sd_post_digits = 1, sd_change_digits = 1)
  expect_true(hi$recalculated_lower > 1)
  # change SD above sd_pre + sd_post -> r < -1 throughout
  lo <- recalc_prepost_r(8, 5, 13.5,
                         sd_pre_digits = 1, sd_post_digits = 1, sd_change_digits = 1)
  expect_true(lo$recalculated_upper < -1)
  # a feasible value -> interval inside [-1, 1]
  ok <- recalc_prepost_r(5, 8, 6,
                         sd_pre_digits = 1, sd_post_digits = 1, sd_change_digits = 1)
  expect_true(ok$recalculated_lower >= -1 && ok$recalculated_upper <= 1)
})

test_that("recalc_prepost_r interval widens with coarser reporting", {
  fine   <- recalc_prepost_r(6.77, 7.45, 13.01, 4, 4, 4)
  coarse <- recalc_prepost_r(6.8,  7.5,  13.0,  1, 1, 1)
  w <- function(x) x$recalculated_upper - x$recalculated_lower
  expect_gt(w(coarse), w(fine))
  expect_gte(w(fine), 0)
})

test_that("recalc_prepost_r consistency check against a reported r", {
  inside  <- recalc_prepost_r(6.77, 7.45, 13.01, 2, 2, 2, r = -0.67, r_digits = 2)
  outside <- recalc_prepost_r(6.77, 7.45, 13.01, 2, 2, 2, r =  0.50, r_digits = 2)
  expect_true(inside$consistent)
  expect_false(outside$consistent)
})

test_that("recalc_prepost_r requires digits", {
  expect_error(recalc_prepost_r(6.77, 7.45, 13.01), "digits")
})

test_that("recalc_change_sd_from_r inverts recalc_prepost_r", {
  # text SDs + r 0.458 should reproduce change SD ~ 13.01
  res <- recalc_change_sd_from_r(14.58, 7.79, 0.458,
                                 sd_pre_digits = 2, sd_post_digits = 2, r_digits = 3)
  mid <- (res$recalculated_lower + res$recalculated_upper) / 2
  expect_equal(mid, 13.01, tolerance = 0.02)
})

test_that("recalc_prepost_r_from_f reproduces Pu et al. (2026)", {
  res <- recalc_prepost_r_from_f(
    f = 117.055,
    m1b = 24.71, sd1b = 4.137, m1p = 8.96,  sd1p = 5.237, n1 = 52,
    m2b = 24.81, sd2b = 3.774, m2p = 15.95, sd2p = 5.714, n2 = 26,
    f_digits = 3, m_digits = 2, sd_digits = 3)
  mid <- (res$recalculated_lower + res$recalculated_upper) / 2
  expect_equal(mid, 0.885, tolerance = 0.01)
  expect_true(res$recalculated_lower > 0 && res$recalculated_upper < 1)
})

test_that("recalc_prepost_r_from_f flags an F incompatible with the SDs", {
  res <- recalc_prepost_r_from_f(
    f = 10000,
    m1b = 24.71, sd1b = 4.137, m1p = 8.96,  sd1p = 5.237, n1 = 52,
    m2b = 24.81, sd2b = 3.774, m2p = 15.95, sd2p = 5.714, n2 = 26,
    f_digits = 1, m_digits = 2, sd_digits = 3)
  expect_true(res$recalculated_lower > 1)   # implied r > 1 -> impossible
})

test_that("recalc_partial_eta_from_f reproduces F <-> partial eta^2", {
  res <- recalc_partial_eta_from_f(117.055, df_effect = 1, df_error = 76,
                                   f_digits = 3, eta = 0.606, eta_digits = 3)
  mid <- (res$recalculated_lower + res$recalculated_upper) / 2
  expect_equal(mid, 0.606, tolerance = 1e-3)
  expect_true(res$consistent)
  # several reported values reconcile
  Fv  <- c(32.713, 51.214, 25.263)
  eta <- c(0.301,  0.403,  0.249)
  for (i in seq_along(Fv)) {
    r <- recalc_partial_eta_from_f(Fv[i], 1, 76, f_digits = 3, eta = eta[i], eta_digits = 3)
    expect_true(r$consistent)
  }
})
