#' Recalculate F from R^2, N, k
#'
#' Identity: \eqn{F = (R^2/k) / ((1-R^2)/(N-k-1))}. Composes with
#' \code{recalc_p_from_f()} to recover model p. Most common chain: paper omits
#' F but reports R^2, N, k.
#'
#' @param r2,n,k Numeric.
#' @param f Numeric. Reported F (optional).
#' @param r2_digits,n_digits,k_digits,f_digits Integer.
#' @return One-row tibble.
#' @examples
#' recalc_f_from_r2(r2 = 0.41, n = 1923, k = 2, f = 665.4)
#' @export
#' @param rounding See \code{\link{recalc_rounding}} for the accepted values.
recalc_f_from_r2 <- function(
  r2,
  n,
  k,
  f = NULL,
  r2_digits = 2,
  n_digits = 0,
  k_digits = 0,
  f_digits = 2,
  rounding = "either"
) {
  recomp <- propagate_intervals(
    fn = function(r2, n, k) (r2 / k) / ((1 - r2) / (n - k - 1)),
    inputs = list(
      r2 = interval_from_digits(r2, r2_digits, rounding = rounding),
      n = interval_from_digits(n, n_digits, rounding = rounding),
      k = interval_from_digits(k, k_digits, rounding = rounding)
    )
  )
  recalc_result(
    "B1: F = (R^2/k) / ((1-R^2)/(N-k-1))",
    f,
    reported_interval(f, f_digits, rounding = rounding),
    recomp
  )
}

#' Recalculate R^2 from F, N, k
#'
#' Identity: \eqn{R^2 = Fk / (Fk + N - k - 1)}.
#'
#' @param f,n,k Numeric.
#' @param r2 Numeric. Reported R^2 (optional).
#' @param f_digits,n_digits,k_digits,r2_digits Integer.
#' @return One-row tibble.
#' @examples
#' recalc_r2_from_f(f = 665.4, n = 1923, k = 2, r2 = 0.41)
#' @export
#' @param rounding See \code{\link{recalc_rounding}} for the accepted values.
recalc_r2_from_f <- function(
  f,
  n,
  k,
  r2 = NULL,
  f_digits = 2,
  n_digits = 0,
  k_digits = 0,
  r2_digits = 2,
  rounding = "either"
) {
  recomp <- propagate_intervals(
    fn = function(f, n, k) (f * k) / (f * k + (n - k - 1)),
    inputs = list(
      f = interval_from_digits(f, f_digits, rounding = rounding),
      n = interval_from_digits(n, n_digits, rounding = rounding),
      k = interval_from_digits(k, k_digits, rounding = rounding)
    )
  )
  recalc_result(
    "B1: R^2 = Fk / (Fk + N-k-1)",
    r2,
    reported_interval(r2, r2_digits, rounding = rounding),
    recomp
  )
}

#' Recalculate omnibus model p from F and df
#'
#' Identity: \eqn{p_\text{model} = 1 - F_F(F; k, N-k-1)}. Composes with
#' \code{recalc_f_from_r2()}: when F is not reported, recompute it from R^2
#' first, then call this function. See also \code{recalc_p_model_from_r2()} for
#' the chained version with rounding propagated across both stages.
#'
#' @param f,df1,df2 Numeric.
#' @param p Numeric. Reported model p (optional).
#' @param p_op One of \code{"eq"}, \code{"lt"}, \code{"le"}, \code{"gt"},
#'   \code{"ge"}.
#' @param f_digits,df1_digits,df2_digits,p_digits Integer.
#' @return One-row tibble.
#' @examples
#' recalc_p_from_f(f = 665.4, df1 = 2, df2 = 1920, p = 0.001, p_op = "lt")
#' @export
#' @param rounding See \code{\link{recalc_rounding}} for the accepted values.
recalc_p_from_f <- function(
  f,
  df1,
  df2,
  p = NULL,
  p_op = "eq",
  f_digits = 2,
  df1_digits = 0,
  df2_digits = 0,
  p_digits = 3,
  rounding = "either"
) {
  recomp <- propagate_intervals(
    fn = function(f, df1, df2) {
      stats::pf(f, df1 = df1, df2 = df2, lower.tail = FALSE)
    },
    inputs = list(
      f = interval_from_digits(f, f_digits, rounding = rounding),
      df1 = interval_from_digits(df1, df1_digits, rounding = rounding),
      df2 = interval_from_digits(df2, df2_digits, rounding = rounding)
    )
  )
  recalc_result(
    "B2: p_model = 1 - F_F(F; df1, df2)",
    p,
    reported_interval(p, p_digits, op = p_op, rounding = rounding),
    recomp
  )
}

#' Recalculate adjusted R^2 from R^2, N, k
#'
#' Identity: \eqn{R^2_\text{adj} = 1 - (1 - R^2)(N-1)/(N-k-1)}. Also the
#' diagnostic for R^2 / adjusted R^2 label confusion; see
#' \code{diagnose_r2_label()}.
#'
#' @param r2,n,k Numeric.
#' @param adj_r2 Numeric. Reported adjusted R^2 (optional).
#' @param r2_digits,n_digits,k_digits,adj_r2_digits Integer.
#' @return One-row tibble.
#' @examples
#' recalc_adj_r2(r2 = 0.41, n = 1923, k = 2, adj_r2 = 0.41)
#' @export
#' @param rounding See \code{\link{recalc_rounding}} for the accepted values.
recalc_adj_r2 <- function(
  r2,
  n,
  k,
  adj_r2 = NULL,
  r2_digits = 2,
  n_digits = 0,
  k_digits = 0,
  adj_r2_digits = 2,
  rounding = "either"
) {
  recomp <- propagate_intervals(
    fn = function(r2, n, k) 1 - (1 - r2) * (n - 1) / (n - k - 1),
    inputs = list(
      r2 = interval_from_digits(r2, r2_digits, rounding = rounding),
      n = interval_from_digits(n, n_digits, rounding = rounding),
      k = interval_from_digits(k, k_digits, rounding = rounding)
    )
  )
  recalc_result(
    "B3: R^2_adj = 1 - (1-R^2)(N-1)/(N-k-1)",
    adj_r2,
    reported_interval(adj_r2, adj_r2_digits, rounding = rounding),
    recomp
  )
}

#' Recalculate R^2 from standardized betas and zero-order correlations
#'
#' Identity: \eqn{R^2 = \sum_i \beta_i r_{Yi}}. Inputs are vectors; rounding
#' is propagated by enumerating \eqn{2^{2k}} corner combinations.
#'
#' @param betas Numeric vector. Standardized regression coefficients.
#' @param r_y Numeric vector. Zero-order correlations with Y (same length as
#'   \code{betas}).
#' @param r2 Numeric. Reported model R^2 (optional).
#' @param betas_digits,r_y_digits,r2_digits Integer.
#' @return One-row tibble.
#' @examples
#' recalc_r2_from_betas_corrs(betas = c(-0.39, 0.47),
#'                            r_y = c(-0.53, 0.61), r2 = 0.41)
#' @export
#' @param rounding See \code{\link{recalc_rounding}} for the accepted values.
recalc_r2_from_betas_corrs <- function(
  betas,
  r_y,
  r2 = NULL,
  betas_digits = 2,
  r_y_digits = 2,
  r2_digits = 2,
  rounding = "either"
) {
  stopifnot(length(betas) == length(r_y))
  k <- length(betas)
  inputs <- c(
    setNames(
      lapply(betas, interval_from_digits, betas_digits, rounding = rounding),
      paste0("b", seq_len(k))
    ),
    setNames(
      lapply(r_y, interval_from_digits, r_y_digits, rounding = rounding),
      paste0("r", seq_len(k))
    )
  )
  fn <- function(...) {
    args <- list(...)
    bs <- unlist(args[paste0("b", seq_len(k))])
    rs <- unlist(args[paste0("r", seq_len(k))])
    sum(bs * rs)
  }
  recomp <- propagate_intervals(fn, inputs)
  recalc_result(
    "B4: R^2 = sum_i beta_i * r_Yi",
    r2,
    reported_interval(r2, r2_digits, rounding = rounding),
    recomp
  )
}

#' Recalculate R^2 from a sequence of hierarchical Delta R^2 blocks
#'
#' Identity: \eqn{R^2_\text{full} = \sum_j \Delta R^2_j}. Holds by construction
#' in any hierarchical regression. Composes with \code{recalc_semipartial_r2_from_t()}:
#' each per-step Delta R^2 should also equal sr_i^2 for the predictor added that
#' step (when the step adds one predictor).
#'
#' @param delta_r2 Numeric vector. Delta R^2 for each block.
#' @param r2 Numeric. Reported full-model R^2 (optional).
#' @param delta_r2_digits,r2_digits Integer.
#' @return One-row tibble.
#' @examples
#' recalc_r2_from_blocks(delta_r2 = c(0.18, 0.23), r2 = 0.41)
#' @export
#' @param rounding See \code{\link{recalc_rounding}} for the accepted values.
recalc_r2_from_blocks <- function(
  delta_r2,
  r2 = NULL,
  delta_r2_digits = 2,
  r2_digits = 2,
  rounding = "either"
) {
  k <- length(delta_r2)
  inputs <- setNames(
    lapply(
      delta_r2,
      interval_from_digits,
      delta_r2_digits,
      rounding = rounding
    ),
    paste0("d", seq_len(k))
  )
  fn <- function(...) sum(unlist(list(...)))
  recomp <- propagate_intervals(fn, inputs)
  recalc_result(
    "B5: R^2 = sum_j Delta R^2_j",
    r2,
    reported_interval(r2, r2_digits, rounding = rounding),
    recomp
  )
}
