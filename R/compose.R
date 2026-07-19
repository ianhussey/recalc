#' Recover SE from a Wald CI and df
#'
#' Inverse of A3: \eqn{\text{SE} = (\text{CI}_U - \text{CI}_L)/(2 t_{\text{crit}})}.
#' Useful when the paper reports b and CI but not SE; the recovered SE feeds
#' \code{recalc_t_from_b_se()} (to get t) and \code{recalc_p_from_t_df()} (to
#' get p).
#'
#' @param ci Length-2 numeric \code{c(lower, upper)}.
#' @param df Numeric.
#' @param level Numeric. Default 0.95.
#' @param ci_digits,df_digits Integer.
#' @return One-row tibble. \code{reported} is NA; \code{recalculated_*} is the SE
#'   interval.
#' @examples
#' recalc_se_from_ci(ci = c(-6.49, -4.20), df = 30)
#' @export
#' @param rounding See \code{\link{recalc_rounding}} for the accepted values.
recalc_se_from_ci <- function(
  ci,
  df,
  level = 0.95,
  ci_digits = 2,
  df_digits = 0,
  rounding = "either"
) {
  alpha <- 1 - level
  recomp <- propagate_intervals(
    fn = function(lo, hi, df) (hi - lo) / (2 * stats::qt(1 - alpha / 2, df)),
    inputs = list(
      lo = interval_from_digits(ci[1], ci_digits, rounding = rounding),
      hi = interval_from_digits(ci[2], ci_digits, rounding = rounding),
      df = interval_from_digits(df, df_digits, rounding = rounding)
    )
  )
  recalc_result(
    "D: SE = (CI_U - CI_L) / (2 t_crit)",
    NA_real_,
    c(NA_real_, NA_real_),
    recomp
  )
}

#' Recover model p from R^2 directly (composition: B1 -> B2)
#'
#' When F is not reported but R^2, N, and k are, model p can be computed by
#' chaining B1 (F from R^2) into B2 (p from F). This helper runs the chain
#' as a single propagation so rounding intervals are handled correctly across
#' both stages.
#'
#' @param r2,n,k Numeric.
#' @param p_model Numeric. Reported model p (optional).
#' @param p_op One of \code{"eq"}, \code{"lt"}, \code{"le"}, \code{"gt"},
#'   \code{"ge"}.
#' @param r2_digits,n_digits,k_digits,p_model_digits Integer.
#' @return One-row tibble.
#' @examples
#' recalc_p_model_from_r2(r2 = 0.41, n = 1923, k = 2,
#'                        p_model = 0.001, p_op = "lt")
#' @export
#' @param rounding See \code{\link{recalc_rounding}} for the accepted values.
recalc_p_model_from_r2 <- function(
  r2,
  n,
  k,
  p_model = NULL,
  p_op = "eq",
  r2_digits = 2,
  n_digits = 0,
  k_digits = 0,
  p_model_digits = 3,
  rounding = "either"
) {
  recomp <- propagate_intervals(
    fn = function(r2, n, k) {
      f <- (r2 / k) / ((1 - r2) / (n - k - 1))
      stats::pf(f, df1 = k, df2 = n - k - 1, lower.tail = FALSE)
    },
    inputs = list(
      r2 = interval_from_digits(r2, r2_digits, rounding = rounding),
      n = interval_from_digits(n, n_digits, rounding = rounding),
      k = interval_from_digits(k, k_digits, rounding = rounding)
    )
  )
  recalc_result(
    "F: p_model from R^2 (B1 -> B2)",
    p_model,
    reported_interval(p_model, p_model_digits, op = p_op, rounding = rounding),
    recomp
  )
}
