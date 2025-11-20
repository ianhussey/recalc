#' Calculate correlation coefficients corresponding to significance levels
#'
#' This function calculates the correlation coefficients that correspond to specific
#' significance levels (alpha levels) for a given sample size. This can be used to see
#' if correlation stars (or sample sizes) are consistent across a correlation table.
#' It uses Fisher's z transformation to convert alpha levels into correlation coefficients.
#'
#' @param n An integer representing the sample size.
#' @param alpha_levels A numeric vector of significance levels (default is c(0.05, 0.01, 0.001)).
#' @return A tibble with two columns: `p_values`, which indicates the significance levels,
#' and `r_values`, which are the corresponding correlation coefficients.
#' @examples
#' get_correlation_thresholds(30)
#' get_correlation_thresholds(100, alpha_levels = c(0.05, 0.01))
#' @export

get_correlation_thresholds <- function(n, alpha_levels = c(0.05, 0.01, 0.001)) {
  # Calculate z scores for the given alpha levels (two-tailed)
  z_scores <- stats::qnorm(1 - alpha_levels / 2)

  # Convert z scores to Fisher z' values
  fisher_z_values <- z_scores * sqrt(1 / (n - 3))

  # Convert Fisher z' values back to correlation coefficients
  r_values <- tanh(fisher_z_values)

  tibble::tibble(p_values = paste("p <=", alpha_levels), r_values = r_values)
}
