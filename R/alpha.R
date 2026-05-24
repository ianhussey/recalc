#' Recalculate Cronbach's alpha from a correlation matrix and SD vector
#'
#' Identity: \eqn{\alpha = \frac{p}{p-1}\left(1 - \sum_j s_j^2 / V\right)} where
#' \eqn{V = \sum_j s_j^2 + 2 \sum_{j<k} r_{jk} s_j s_k} is the total variance of
#' the summed score implied by the reported item correlation matrix
#' \code{cor_mat} and SD vector \code{sd_vector}. Rounding is propagated by
#' enumerating \eqn{2^{p(p+1)/2}} corner combinations (\eqn{p} item SDs plus
#' \eqn{p(p-1)/2} unique correlations); cost grows quickly past about
#' \eqn{p = 6} items.
#'
#' Mismatch with the reported alpha is a forensic signal: \code{cor_mat} and
#' \code{sd_vector} were not computed from the same data, or \code{alpha} was
#' computed from a different covariance / scoring rule (e.g. with imputed
#' items).
#'
#' @param cor_mat Numeric. \eqn{p \times p} item correlation matrix.
#' @param sd_vector Numeric. Length-\eqn{p} vector of item SDs.
#' @param alpha Numeric. Reported Cronbach's alpha (optional).
#' @param cor_mat_digits,sd_vector_digits,alpha_digits Integer. Number of
#'   decimal places each value was reported to. No defaults: the caller must
#'   specify the precision of every value they supply. \code{alpha_digits} is
#'   required only when \code{alpha} is provided.
#' @return One-row tibble.
#' @examples
#' R <- matrix(c(1.00, 0.50, 0.40,
#'               0.50, 1.00, 0.60,
#'               0.40, 0.60, 1.00), 3, 3)
#' recalc_alpha_from_cor_sd(cor_mat = R, sd_vector = c(1.20, 1.10, 0.90),
#'                          alpha = 0.78,
#'                          cor_mat_digits = 2, sd_vector_digits = 2,
#'                          alpha_digits = 2)
#' @export
recalc_alpha_from_cor_sd <- function(cor_mat, sd_vector, alpha = NULL,
                                     cor_mat_digits = NULL,
                                     sd_vector_digits = NULL,
                                     alpha_digits = NULL, rounding = "either") {
  require_digits(cor_mat_digits = cor_mat_digits,
                 sd_vector_digits = sd_vector_digits)
  if (!is.null(alpha)) require_digits(alpha_digits = alpha_digits)

  stopifnot(is.matrix(cor_mat),
            nrow(cor_mat) == ncol(cor_mat),
            length(sd_vector) == nrow(cor_mat))
  p <- length(sd_vector)
  if (p < 2L) stop("alpha requires at least 2 items")

  s_names <- paste0("s", seq_len(p))
  idx <- which(upper.tri(cor_mat), arr.ind = TRUE)
  r_names <- sprintf("r_%d_%d", idx[, 1], idx[, 2])

  inputs <- c(
    setNames(lapply(sd_vector, interval_from_digits, sd_vector_digits, rounding = rounding),
             s_names),
    setNames(lapply(seq_len(nrow(idx)),
                    function(i) interval_from_digits(
                      cor_mat[idx[i, 1], idx[i, 2]], cor_mat_digits,
                      rounding = rounding)),
             r_names)
  )

  fn <- function(...) {
    args <- list(...)
    ss <- unlist(args[s_names])
    rr <- unlist(args[r_names])
    sum_ss <- sum(ss^2)
    V <- sum_ss + 2 * sum(rr * ss[idx[, 1]] * ss[idx[, 2]])
    (p / (p - 1)) * (1 - sum_ss / V)
  }
  recomp <- propagate_intervals(fn, inputs)
  recalc_result("D1: alpha = (p/(p-1))(1 - sum s_j^2 / V)",
                alpha, reported_interval(alpha, alpha_digits, rounding = rounding), recomp)
}
