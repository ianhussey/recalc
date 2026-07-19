#' Diagnose whether a coefficient column is standardized or unstandardized
#'
#' Tests two interpretations of a coefficient column reported as "beta":
#' (A) the values are true betas; (B) the values are unstandardized b's
#' mislabeled as beta. Returns the recalculated-R^2 interval under each
#' interpretation so the user can see which (if either) is consistent with
#' the reported R^2.
#'
#' If A is consistent and B is not, the labeling is correct. If B is consistent
#' and A is not, the column is mislabeled. If neither is consistent, the gap
#' is structural.
#'
#' @param coef_reported Numeric vector. Reported values from the column
#'   labeled "beta" (or "b").
#' @param r_y Numeric vector. Zero-order correlations of predictors with Y.
#' @param sd_x Numeric vector. Predictor SDs (one per predictor).
#' @param sd_y Numeric scalar. Outcome SD.
#' @param r2 Numeric. Reported model R^2.
#' @param coef_digits,r_y_digits,sd_x_digits,sd_y_digits,r2_digits Integer.
#' @return Two-row tibble: one row per interpretation.
#' @examples
#' diagnose_beta_label(
#'   coef_reported = c(-0.39, 0.47),
#'   r_y = c(-0.53, 0.61),
#'   sd_x = c(1.2, 0.5),
#'   sd_y = 1.4,
#'   r2 = 0.41,
#'   sd_x_digits = 1, sd_y_digits = 1
#' )
#' @export
#' @param rounding See \code{\link{recalc_rounding}} for the accepted values.
diagnose_beta_label <- function(
  coef_reported,
  r_y,
  sd_x,
  sd_y,
  r2,
  coef_digits = 2,
  r_y_digits = 2,
  sd_x_digits = 1,
  sd_y_digits = 1,
  r2_digits = 2,
  rounding = "either"
) {
  stopifnot(
    length(coef_reported) == length(r_y),
    length(sd_x) == length(coef_reported)
  )
  k <- length(coef_reported)

  # Interpretation A: column is true beta - apply B4 directly
  res_A <- recalc_r2_from_betas_corrs(
    betas = coef_reported,
    r_y = r_y,
    r2 = r2,
    betas_digits = coef_digits,
    r_y_digits = r_y_digits,
    r2_digits = r2_digits
  )
  res_A$check <- "E1 (Interp A): column is true beta"

  # Interpretation B: column is b mislabeled as beta - convert and apply B4
  inputs_B <- c(
    setNames(
      lapply(
        coef_reported,
        interval_from_digits,
        coef_digits,
        rounding = rounding
      ),
      paste0("b", seq_len(k))
    ),
    setNames(
      lapply(sd_x, interval_from_digits, sd_x_digits, rounding = rounding),
      paste0("sx", seq_len(k))
    ),
    setNames(
      lapply(r_y, interval_from_digits, r_y_digits, rounding = rounding),
      paste0("r", seq_len(k))
    ),
    list(sy = interval_from_digits(sd_y, sd_y_digits, rounding = rounding))
  )
  fn_B <- function(...) {
    args <- list(...)
    bs <- unlist(args[paste0("b", seq_len(k))])
    sxs <- unlist(args[paste0("sx", seq_len(k))])
    rs <- unlist(args[paste0("r", seq_len(k))])
    sy <- args[["sy"]]
    sum((bs * sxs / sy) * rs)
  }
  recomp_B <- propagate_intervals(fn_B, inputs_B)
  res_B <- recalc_result(
    "E1 (Interp B): column is b mislabeled as beta",
    r2,
    reported_interval(r2, r2_digits, rounding = rounding),
    recomp_B
  )

  dplyr::bind_rows(res_A, res_B)
}

#' Diagnose whether a reported "R^2" is actually adjusted R^2
#'
#' Tests two interpretations: (A) the value is raw R^2 - apply B4 directly;
#' (B) the value is adjusted R^2 - back out the implied raw R^2 via
#' \eqn{R^2_\text{raw} = 1 - (1 - R^2_\text{adj})(N - k - 1)/(N - 1)} and
#' test that against B4.
#'
#' @param r2_reported Numeric.
#' @param betas Numeric vector.
#' @param r_y Numeric vector.
#' @param n,k Numeric.
#' @param r2_digits,betas_digits,r_y_digits,n_digits,k_digits Integer.
#' @return Two-row tibble.
#' @examples
#' diagnose_r2_label(
#'   r2_reported = 0.41,
#'   betas = c(-0.39, 0.47), r_y = c(-0.53, 0.61),
#'   n = 1923, k = 2
#' )
#' @export
#' @param rounding See \code{\link{recalc_rounding}} for the accepted values.
diagnose_r2_label <- function(
  r2_reported,
  betas,
  r_y,
  n,
  k,
  r2_digits = 2,
  betas_digits = 2,
  r_y_digits = 2,
  n_digits = 0,
  k_digits = 0,
  rounding = "either"
) {
  stopifnot(length(betas) == length(r_y))
  k_pred <- length(betas)

  # B4-implied R^2 interval (independent of which label is correct)
  inputs_br <- c(
    setNames(
      lapply(betas, interval_from_digits, betas_digits, rounding = rounding),
      paste0("b", seq_len(k_pred))
    ),
    setNames(
      lapply(r_y, interval_from_digits, r_y_digits, rounding = rounding),
      paste0("r", seq_len(k_pred))
    )
  )
  fn_br <- function(...) {
    args <- list(...)
    bs <- unlist(args[paste0("b", seq_len(k_pred))])
    rs <- unlist(args[paste0("r", seq_len(k_pred))])
    sum(bs * rs)
  }
  recomp_br <- propagate_intervals(fn_br, inputs_br)

  # Interpretation A: reported value is raw R^2
  res_A <- recalc_result(
    "E2 (Interp A): reported is raw R^2",
    r2_reported,
    reported_interval(r2_reported, r2_digits, rounding = rounding),
    recomp_br
  )

  # Interpretation B: reported value is adjusted R^2 - invert to raw R^2
  raw_int <- propagate_intervals(
    fn = function(adj, n, k) 1 - (1 - adj) * (n - k - 1) / (n - 1),
    inputs = list(
      adj = interval_from_digits(r2_reported, r2_digits, rounding = rounding),
      n = interval_from_digits(n, n_digits, rounding = rounding),
      k = interval_from_digits(k, k_digits, rounding = rounding)
    )
  )
  res_B <- tibble::tibble(
    check = "E2 (Interp B): reported is adjusted R^2 (back out raw)",
    reported = NA_real_,
    reported_lower = raw_int[["lower"]],
    reported_upper = raw_int[["upper"]],
    recalculated_lower = recomp_br[["lower"]],
    recalculated_upper = recomp_br[["upper"]],
    consistent = (recomp_br[["upper"]] >= raw_int[["lower"]]) &
      (recomp_br[["lower"]] <= raw_int[["upper"]])
  )

  dplyr::bind_rows(res_A, res_B)
}

#' Compare one-tailed and two-tailed p interpretations for a reported p
#'
#' A reported p that fails A2 under the default two-tailed assumption may be
#' consistent with a one-tailed test. Runs A2 under both tail conventions and
#' reports both intervals; if exactly one is consistent, the tail convention
#' is the issue.
#'
#' @param t,df Numeric.
#' @param p Numeric. Reported p.
#' @param p_op,t_digits,df_digits,p_digits As in \code{recalc_p_from_t_df()}.
#' @return Two-row tibble.
#' @examples
#' diagnose_p_tails(t = 1.85, df = 200, p = 0.033)
#' @export
#' @param rounding See \code{\link{recalc_rounding}} for the accepted values.
diagnose_p_tails <- function(
  t,
  df,
  p,
  p_op = "eq",
  t_digits = 2,
  df_digits = 0,
  p_digits = 3,
  rounding = "either"
) {
  one <- recalc_p_from_t_df(
    t,
    df,
    p = p,
    p_op = p_op,
    two_tailed = FALSE,
    t_digits = t_digits,
    df_digits = df_digits,
    p_digits = p_digits
  )
  one$check <- "E3 (Interp A): one-tailed p"
  two <- recalc_p_from_t_df(
    t,
    df,
    p = p,
    p_op = p_op,
    two_tailed = TRUE,
    t_digits = t_digits,
    df_digits = df_digits,
    p_digits = p_digits
  )
  two$check <- "E3 (Interp B): two-tailed p"
  dplyr::bind_rows(one, two)
}
