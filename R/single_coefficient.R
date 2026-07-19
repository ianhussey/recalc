#' Recalculate t from a coefficient and its standard error
#'
#' Identity: \eqn{t_i = b_i / \text{SE}_i}. Catches transcription errors and
#' sign flips. Composes with \code{recalc_p_from_t_df()} to recover p.
#'
#' @param b,se Numeric.
#' @param t Numeric. Reported t (optional).
#' @param b_digits,se_digits,t_digits Integer.
#' @return One-row tibble.
#' @examples
#' recalc_t_from_b_se(b = -5.34, se = 0.56, t = -9.56)
#' @export
#' @param rounding See \code{\link{recalc_rounding}} for the accepted values.
recalc_t_from_b_se <- function(
  b,
  se,
  t = NULL,
  b_digits = 2,
  se_digits = 2,
  t_digits = 2,
  rounding = "either"
) {
  recomp <- propagate_intervals(
    fn = function(b, se) b / se,
    inputs = list(
      b = interval_from_digits(b, b_digits, rounding = rounding),
      se = interval_from_digits(se, se_digits, rounding = rounding)
    )
  )
  recalc_result(
    "A1: t = b / SE",
    t,
    reported_interval(t, t_digits, rounding = rounding),
    recomp
  )
}

#' Recalculate a coefficient p from its t and residual df
#'
#' Identity: \eqn{p_i = 2(1 - F_t(|t_i|; \text{df}))}. Use
#' \code{p_op = "lt"} for reports of the form \code{p < .001}. Composes with
#' \code{recalc_t_from_b_se()} to recover t first if only b and SE are reported.
#'
#' @param t,df Numeric.
#' @param p Numeric.
#' @param p_op One of \code{"eq"}, \code{"lt"}, \code{"le"}, \code{"gt"},
#'   \code{"ge"}.
#' @param two_tailed Logical.
#' @param t_digits,df_digits,p_digits Integer.
#' @return One-row tibble.
#' @examples
#' recalc_p_from_t_df(t = -9.56, df = 30, p = 0.001, p_op = "lt")
#' @export
#' @param rounding See \code{\link{recalc_rounding}} for the accepted values.
recalc_p_from_t_df <- function(
  t,
  df,
  p = NULL,
  p_op = "eq",
  two_tailed = TRUE,
  t_digits = 2,
  df_digits = 0,
  p_digits = 3,
  rounding = "either"
) {
  fn <- function(t, df) {
    factor <- if (two_tailed) 2 else 1
    factor * stats::pt(abs(t), df = df, lower.tail = FALSE)
  }
  recomp <- propagate_intervals(
    fn,
    inputs = list(
      t = interval_from_digits(t, t_digits, rounding = rounding),
      df = interval_from_digits(df, df_digits, rounding = rounding)
    )
  )
  recalc_result(
    "A2: p = 2(1 - F_t(|t|, df))",
    p,
    reported_interval(p, p_digits, op = p_op, rounding = rounding),
    recomp
  )
}

#' Recalculate a Wald CI from b, SE, and df
#'
#' Identity: \eqn{\text{CI}_{1-\alpha} = b_i \pm t_{1-\alpha/2, \text{df}} \cdot \text{SE}_i}.
#' Useful in reverse: when SE is missing, \code{recalc_se_from_ci()} recovers
#' it from the CI, which then feeds \code{recalc_t_from_b_se()} and
#' \code{recalc_p_from_t_df()}.
#'
#' @param b,se,df Numeric.
#' @param level Numeric. Default 0.95.
#' @param ci Optional length-2 numeric \code{c(lower, upper)}.
#' @param b_digits,se_digits,df_digits,ci_digits Integer.
#' @return Two-row tibble.
#' @examples
#' recalc_ci_from_b_se(b = -5.34, se = 0.56, df = 30,
#'                     ci = c(-6.49, -4.20))
#' @export
#' @param rounding See \code{\link{recalc_rounding}} for the accepted values.
recalc_ci_from_b_se <- function(
  b,
  se,
  df,
  level = 0.95,
  ci = NULL,
  b_digits = 2,
  se_digits = 2,
  df_digits = 0,
  ci_digits = 2,
  rounding = "either"
) {
  alpha <- 1 - level
  fn_lower <- function(b, se, df) b - stats::qt(1 - alpha / 2, df) * se
  fn_upper <- function(b, se, df) b + stats::qt(1 - alpha / 2, df) * se
  inputs <- list(
    b = interval_from_digits(b, b_digits, rounding = rounding),
    se = interval_from_digits(se, se_digits, rounding = rounding),
    df = interval_from_digits(df, df_digits, rounding = rounding)
  )
  rec_lo <- propagate_intervals(fn_lower, inputs)
  rec_hi <- propagate_intervals(fn_upper, inputs)
  dplyr::bind_rows(
    recalc_result(
      "A3: CI lower = b - t_crit SE",
      if (is.null(ci)) NULL else ci[1],
      reported_interval(
        if (is.null(ci)) NULL else ci[1],
        ci_digits,
        rounding = rounding
      ),
      rec_lo
    ),
    recalc_result(
      "A3: CI upper = b + t_crit SE",
      if (is.null(ci)) NULL else ci[2],
      reported_interval(
        if (is.null(ci)) NULL else ci[2],
        ci_digits,
        rounding = rounding
      ),
      rec_hi
    )
  )
}

#' Recalculate beta from an unstandardized coefficient and SDs
#'
#' Identity: \eqn{\beta_i = b_i \cdot \text{SD}_{X_i} / \text{SD}_Y}. The
#' diagnostic for b/beta label confusion; see \code{diagnose_beta_label()}.
#'
#' @param b,sd_x,sd_y Numeric.
#' @param beta Numeric.
#' @param b_digits,sd_x_digits,sd_y_digits,beta_digits Integer.
#' @return One-row tibble.
#' @examples
#' recalc_beta_from_b(b = -5.34, sd_x = 0.98, sd_y = 6.03, beta = -0.87)
#' @export
#' @param rounding See \code{\link{recalc_rounding}} for the accepted values.
recalc_beta_from_b <- function(
  b,
  sd_x,
  sd_y,
  beta = NULL,
  b_digits = 2,
  sd_x_digits = 2,
  sd_y_digits = 2,
  beta_digits = 2,
  rounding = "either"
) {
  recomp <- propagate_intervals(
    fn = function(b, sd_x, sd_y) b * sd_x / sd_y,
    inputs = list(
      b = interval_from_digits(b, b_digits, rounding = rounding),
      sd_x = interval_from_digits(sd_x, sd_x_digits, rounding = rounding),
      sd_y = interval_from_digits(sd_y, sd_y_digits, rounding = rounding)
    )
  )
  recalc_result(
    "A4: beta = b * SD_X / SD_Y",
    beta,
    reported_interval(beta, beta_digits, rounding = rounding),
    recomp
  )
}

#' Upper bound on |t| from standardized beta, R^2, N, k
#'
#' Identity: \eqn{|t_i| \leq |\beta_i| \sqrt{(N-k-1)/(1-R^2)}}, achieved at
#' \eqn{\text{VIF}_i = 1}. A reported |t| above this bound is incompatible with
#' the reported beta, R^2, N, and k.
#'
#' @param beta,r2,n,k Numeric.
#' @param t Numeric. Reported t (compared in absolute value).
#' @param beta_digits,r2_digits,n_digits,k_digits,t_digits Integer.
#' @return One-row tibble. \code{consistent} is TRUE iff |reported t| does not
#'   exceed the upper bound on |t|.
#' @examples
#' recalc_t_bound_from_beta(beta = 0.30, r2 = 0.20, n = 200, k = 3, t = 4.20)
#' @export
#' @param rounding See \code{\link{recalc_rounding}} for the accepted values.
recalc_t_bound_from_beta <- function(
  beta,
  r2,
  n,
  k,
  t = NULL,
  beta_digits = 2,
  r2_digits = 2,
  n_digits = 0,
  k_digits = 0,
  t_digits = 2,
  rounding = "either"
) {
  recomp <- propagate_intervals(
    fn = function(beta, r2, n, k) abs(beta) * sqrt((n - k - 1) / (1 - r2)),
    inputs = list(
      beta = interval_from_digits(beta, beta_digits, rounding = rounding),
      r2 = interval_from_digits(r2, r2_digits, rounding = rounding),
      n = interval_from_digits(n, n_digits, rounding = rounding),
      k = interval_from_digits(k, k_digits, rounding = rounding)
    )
  )
  out <- recalc_result(
    "A6 bound: |t| <= |beta| sqrt((N-k-1)/(1-R^2))",
    if (is.null(t)) NULL else abs(t),
    reported_interval(
      if (is.null(t)) NULL else abs(t),
      t_digits,
      rounding = rounding
    ),
    recomp
  )
  if (!is.null(t)) {
    out$consistent <- (abs(t) - 0.5 * 10^(-t_digits)) <= recomp[["upper"]]
  }
  out
}

#' Recalculate squared semipartial (= Delta R^2) from t, R^2, df
#'
#' Identity: \eqn{\text{sr}_i^2 = t_i^2 (1 - R^2) / \text{df}}. The squared
#' semipartial equals Delta R^2 for adding X_i to the model. Useful in
#' hierarchical regressions that report per-step Delta R^2.
#'
#' @param t,r2,df Numeric.
#' @param sr2 Numeric. Reported squared semipartial (optional).
#' @param t_digits,r2_digits,df_digits,sr2_digits Integer.
#' @return One-row tibble.
#' @examples
#' recalc_semipartial_r2_from_t(t = 3.50, r2 = 0.40, df = 200, sr2 = 0.037)
#' @export
#' @param rounding See \code{\link{recalc_rounding}} for the accepted values.
recalc_semipartial_r2_from_t <- function(
  t,
  r2,
  df,
  sr2 = NULL,
  t_digits = 2,
  r2_digits = 2,
  df_digits = 0,
  sr2_digits = 3,
  rounding = "either"
) {
  recomp <- propagate_intervals(
    fn = function(t, r2, df) t^2 * (1 - r2) / df,
    inputs = list(
      t = interval_from_digits(t, t_digits, rounding = rounding),
      r2 = interval_from_digits(r2, r2_digits, rounding = rounding),
      df = interval_from_digits(df, df_digits, rounding = rounding)
    )
  )
  recalc_result(
    "A8: sr_i^2 = t_i^2 (1 - R^2) / df",
    sr2,
    reported_interval(sr2, sr2_digits, rounding = rounding),
    recomp
  )
}
