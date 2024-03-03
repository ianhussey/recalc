test_that("Correct correlation coefficients are returned for n = 30 in a tibble", {
  expected_values <- tibble::tibble(
    p_values = c("p <= 0.05", "p <= 0.01", "p <= 0.001"),
    r_values = c(0.36, 0.46, 0.56)
  )

  result <- get_correlation_thresholds(30)

  expect_equal(round(result$r_values[1], 2), expected_values$r_values[1])
  expect_equal(round(result$r_values[2], 2), expected_values$r_values[2])
  expect_equal(round(result$r_values[3], 2), expected_values$r_values[3])
  expect_equal(result$p_values, expected_values$p_values)
})
