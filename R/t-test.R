#' Recalculate t-test Statistics
#'
#' This function recalculates t-test statistics based on sample means, standard deviations, and sizes. It supports both pooled and Welch's t-test calculations.
#'
#' @param mean1 Mean of the first group.
#' @param sd1 Standard deviation of the first group.
#' @param n1 Sample size of the first group.
#' @param mean2 Mean of the second group.
#' @param sd2 Standard deviation of the second group (if reported separately)
#' @param n2 Sample size of the second group.
#' @param welch Logical; if `TRUE`, uses Welch's t-test, otherwise uses pooled variance t-test.
#' @param round A placeholder for future implementation to round the results. Currently, not implemented.
#' @param alternative A placeholder for future implementation to specify the alternative hypothesis. Currently, not implemented.
#' @param paired Logical; indicates if the t-test should be paired. Currently, not supported.
#' @param round A placeholder for future implementation to round the results. Currently, not implemented.
#' @param alternative A placeholder for future implementation to specify the alternative hypothesis. Currently, not implemented.
#' @return A dataframe containing the t-statistic, degrees of freedom, and p-value.
#' @export
#'
#' @examples
#' recalc_t_test(5, 2, 30, 5.5, 2.5, 30, FALSE)

recalc_t_test <- function(mean1, sd1, n1, mean2, sd2 = NULL, n2, welch = FALSE, paired = FALSE, round = "not implemented", alternative = "not implemented") {

  if (is.null(sd2)) {
    if (welch) "Warning: Welch's t-test requires standard deviation for both groups, returning pooled t-test"
    welch <- FALSE
    sd2 <- sd1
  }

  warning("Currently incorrect - ensure tests work before trusting any of this!")
  if (paired) stop("Paired t-tests not yet supported")
  if (!welch) {
    pooled_sd <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
    t_stat <- (mean1 - mean2) / (pooled_sd * sqrt(1/n1 + 1/n2))
    df <- n1 + n2 - 2
  } else {
    s_delta <- sqrt(sd1^2 / n1 + sd2^2 / n2)
    t_stat <- (mean1 - mean2) / s_delta
    df <- (sd1^2 / n1 + sd2^2 / n2)^2 / ((sd1^2 / n1)^2 / (n1 - 1) + (sd2^2 / n2)^2 / (n2 - 1))
  }
  p_value <- 2 * pt(-abs(t_stat), df)
  return(data.frame(t_statistic = t_stat, degrees_of_freedom = df, p_value = p_value))
}

#' Check Consistency of t-test Results
#'
#' This function checks the consistency of user-provided t-test results against recalculated values.
#' It uses the `recalc_t_test` function to recalculate the statistics and then compares them.
#'
#' @inheritParams recalc_t_test
#' @param t The t-statistic provided by the user.
#' @param df The degrees of freedom provided by the user.
#' @param p The p-value provided by the user.
#' @param welch Logical; if `TRUE`, uses Welch's t-test, otherwise uses pooled variance t-test.
#' @return A dataframe with recalculated statistics and consistency checks for t, degrees of freedom, and p-value.
#' @export
#'
#' @examples
#' check_t_test(5, 2, 30, 5.5, 2.5, 30, 2, 58, 0.05)
check_t_test <- function(mean1, sd1, n1, mean2, sd2 = NULL, n2, t, df, p, welch = FALSE, round = "not implemented", alternative = "not implemented", paired = FALSE) {
  res <- recalc_t_test(mean1, sd1, n1, mean2, sd2, n2, welch = welch)
  res$t_consistent <- with(res, t_statistic < t & t_statistic > t)
  res$df_consistent <- with(res, degrees_of_freedom < df & degrees_of_freedom > df)
  res$p_consistent <- with(res, p_value < p & p_value > p)
  return(res)
}

# Independent t test from summary stats (M, SD, N)
# TODO: like anova_from_summary_stats, this function should accomodate the potential rounding of the reported summary stats
# and return reported, min and max values for t, p, cohen's d, and its CIs.
# TODO: integrate effect_sizes_from_summary_statistics and this function, ie one function that produces t test output (df, p values)
# and also multiple effect sizes (d, g), and also min max values due to rounding

independent_ttest_from_summary_stats <- function(m1, sd1, n1, m2, sd2, n2, sig_level = 0.05){

  # Calculate the t-test statistic
  t_value <- (m1 - m2) / sqrt((sd1^2/n1) + (sd2^2/n2))

  # Calculate the degrees of freedom using the approximation
  df <- ((sd1^2/n1 + sd2^2/n2)^2) / ((sd1^2/n1)^2/(n1-1) + (sd2^2/n2)^2/(n2-1))

  # Obtain p-value for two-tailed test
  p_value <- 2 * (1 - pt(abs(t_value), df))

  # Calculate Cohen's d
  pooled_sd <- sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2) / (n1 + n2 - 2))
  d <- (m1 - m2) / pooled_sd

  # Calculate the standard error of d
  se_d <- sqrt(((n1 + n2) / (n1 * n2)) + (d^2 / (2 * (n1 + n2))))

  # Calculate the t critical value for 95% CI
  t_critical <- qt(1 - (sig_level/2), df = n1 + n2 - 2)

  # Calculate the 95% CI for Cohen's d
  d_ci_lower <- d - t_critical * se_d
  d_ci_upper <- d + t_critical * se_d

  res <-
    data.frame(t = t_value,
               df = df,
               p = p_value,
               cohens_d = d,
               cohens_d_ci_lower =  d_ci_lower,
               cohens_d_ci_upper =  d_ci_upper)

  return(res)
}

#' # Example
#' independent_ttest_from_summary_stats(m1 = 10,
#'                                      sd1 = 2,
#'                                      n1 = 30,
#'                                      m2 = 12,
#'                                      sd2 = 3,
#'                                      n2 = 30)



