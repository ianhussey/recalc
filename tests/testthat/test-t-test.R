# Correct values
tt <- t.test(mpg ~ carb == 4, mtcars, var.equal = TRUE)
tt_welch <- t.test(mpg ~ carb == 4, mtcars, var.equal = FALSE)
t <- round(tt$statistic, 2)
p <- round(tt$p.value, 3)
df <- round(tt$parameter, 2)
t_welch <- round(tt_welch$statistic, 2)
p_welch <- round(tt_welch$p.value, 3)
df_welch <- round(tt_welch$parameter, 2)
m2 <- round(mean(mtcars$mpg[mtcars$carb == 4]), 2)
m1 <- round(mean(mtcars$mpg[mtcars$carb != 4]), 2)
sd2 <- round(sd(mtcars$mpg[mtcars$carb == 4]), 2)
sd1 <- round(sd(mtcars$mpg[mtcars$carb != 4]), 2)

# Test 1: Correctness of recalc_t_test for pooled variance
test_that("recalc_t_test returns correct values for pooled variance", {
  expected <- data.frame(t_statistic = t, degrees_of_freedom = df, p_value = p)
  result <- recalc_t_test(m1, sd1, sum(mtcars$carb == 4), m2, sd2, sum(mtcars$carb != 4), FALSE)
  expect_equal(result$t_statistic, expected$t_statistic, tolerance = 0.01)
  expect_equal(result$degrees_of_freedom, expected$degrees_of_freedom, tolerance = 0.01)
  expect_equal(result$p_value, expected$p_value, tolerance = 0.001)
})

# Test 2: Correctness of recalc_t_test for Welch's t-test
test_that("recalc_t_test returns correct values for Welch's t-test", {
  expected_welch <- data.frame(t_statistic = t_welch, degrees_of_freedom = df_welch, p_value = p_welch)
  result_welch <- recalc_t_test(m1, sd1, sum(mtcars$carb == 4), m2, sd2, sum(mtcars$carb != 4), TRUE)
  expect_equal(result_welch$t_statistic, expected_welch$t_statistic, tolerance = 0.01)
  expect_equal(result_welch$degrees_of_freedom, expected_welch$degrees_of_freedom, tolerance = 0.01)
  expect_equal(result_welch$p_value, expected_welch$p_value, tolerance = 0.001)
})

# Test 3: Error handling for unsupported paired tests
test_that("recalc_t_test correctly handles unsupported paired tests", {
  expect_error(recalc_t_test(m1, sd1, sum(mtcars$carb == 4), m2, sd2, sum(mtcars$carb != 4), FALSE, TRUE))
})

# Test 4: Check_t_test identifies incorrect inputs
test_that("check_t_test identifies inconsistencies", {
  # Simulate incorrect values
  incorrect_t <- t * 1.1
  incorrect_df <- df + 5
  incorrect_p <- p * 1.1

  result <- check_t_test(m1, sd1, sum(mtcars$carb == 4), m2, sd2, sum(mtcars$carb != 4), incorrect_t, incorrect_df, incorrect_p, FALSE)

  # Expect the consistency checks to fail
  expect_false(result$t_consistent)
  expect_false(result$df_consistent)
  expect_false(result$p_consistent)
})
