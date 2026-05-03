#' Two derivations of r_12 from k = 2 standardized normal equations
#'
#' For an internally consistent OLS fit the two derivations of r_12 must agree:
#' \eqn{r_{12} = (r_{Y1} - \beta_1)/\beta_2 = (r_{Y2} - \beta_2)/\beta_1}.
#' Composes with \code{recalc_betas_from_corrs_k2()}: if Table 1 reports r_12,
#' the same value should also satisfy B4 via C2's closed form.
#'
#' @param beta1,beta2 Numeric. Standardized regression coefficients.
#' @param r_y1,r_y2 Numeric. Zero-order correlations with Y.
#' @param beta1_digits,beta2_digits,r_y1_digits,r_y2_digits Integer.
#' @return Two-row tibble: one row per derivation of r_12. The second row uses
#'   the first derivation's interval as its reference to show whether both
#'   derivations agree.
#' @examples
#' recalc_r12_from_normal_eqs(beta1 = -0.39, beta2 = 0.47,
#'                            r_y1 = -0.53, r_y2 = 0.61)
#' @export
recalc_r12_from_normal_eqs <- function(beta1, beta2, r_y1, r_y2,
                                       beta1_digits = 2, beta2_digits = 2,
                                       r_y1_digits = 2, r_y2_digits = 2) {
  inputs <- list(beta1 = interval_from_digits(beta1, beta1_digits),
                 beta2 = interval_from_digits(beta2, beta2_digits),
                 r_y1 = interval_from_digits(r_y1, r_y1_digits),
                 r_y2 = interval_from_digits(r_y2, r_y2_digits))
  recomp1 <- propagate_intervals(
    fn = function(beta1, beta2, r_y1, r_y2) (r_y1 - beta1) / beta2,
    inputs = inputs
  )
  recomp2 <- propagate_intervals(
    fn = function(beta1, beta2, r_y1, r_y2) (r_y2 - beta2) / beta1,
    inputs = inputs
  )
  dplyr::bind_rows(
    recalc_result("C1: r_12 = (r_Y1 - beta_1)/beta_2",
                  NULL, c(NA_real_, NA_real_), recomp1),
    recalc_result("C1: r_12 = (r_Y2 - beta_2)/beta_1",
                  NA_real_,
                  c(recomp1[["lower"]], recomp1[["upper"]]),
                  recomp2)
  )
}

#' Recalculate k = 2 standardized betas from the three correlations
#'
#' Closed-form solution for a k = 2 OLS regression:
#' \eqn{\beta_1 = (r_{Y1} - r_{12} r_{Y2}) / (1 - r_{12}^2)},
#' \eqn{\beta_2 = (r_{Y2} - r_{12} r_{Y1}) / (1 - r_{12}^2)}.
#' Hard cross-check between Table 1 (correlation matrix) and Table 2
#' (regression table) for any k = 2 paper.
#'
#' @param r_y1,r_y2 Numeric. Zero-order correlations with Y.
#' @param r_12 Numeric. Predictor inter-correlation.
#' @param betas Optional length-2 numeric. Reported standardized betas.
#' @param r_y1_digits,r_y2_digits,r_12_digits,betas_digits Integer.
#' @return Two-row tibble: one row per beta.
#' @examples
#' recalc_betas_from_corrs_k2(r_y1 = -0.53, r_y2 = 0.61, r_12 = -0.45,
#'                            betas = c(-0.39, 0.47))
#' @export
recalc_betas_from_corrs_k2 <- function(r_y1, r_y2, r_12, betas = NULL,
                                       r_y1_digits = 2, r_y2_digits = 2,
                                       r_12_digits = 2, betas_digits = 2) {
  inputs <- list(r_y1 = interval_from_digits(r_y1, r_y1_digits),
                 r_y2 = interval_from_digits(r_y2, r_y2_digits),
                 r_12 = interval_from_digits(r_12, r_12_digits))
  rec1 <- propagate_intervals(
    fn = function(r_y1, r_y2, r_12) (r_y1 - r_12 * r_y2) / (1 - r_12^2),
    inputs = inputs
  )
  rec2 <- propagate_intervals(
    fn = function(r_y1, r_y2, r_12) (r_y2 - r_12 * r_y1) / (1 - r_12^2),
    inputs = inputs
  )
  dplyr::bind_rows(
    recalc_result("C2: beta_1 = (r_Y1 - r_12 r_Y2)/(1 - r_12^2)",
                  if (is.null(betas)) NULL else betas[1],
                  reported_interval(if (is.null(betas)) NULL else betas[1],
                                    betas_digits),
                  rec1),
    recalc_result("C2: beta_2 = (r_Y2 - r_12 r_Y1)/(1 - r_12^2)",
                  if (is.null(betas)) NULL else betas[2],
                  reported_interval(if (is.null(betas)) NULL else betas[2],
                                    betas_digits),
                  rec2)
  )
}

#' Lower bound on R^2 from zero-order correlations
#'
#' Identity: \eqn{R^2_\text{full} \geq \max_i r_{Yi}^2}. A reported full-model
#' R^2 below this floor is impossible. Cheap, one-sided, and applicable
#' whenever a correlation table reports correlations of predictors with Y.
#'
#' @param r_y Numeric vector. Zero-order correlations with Y.
#' @param r2 Numeric. Reported model R^2 (optional).
#' @param r_y_digits,r2_digits Integer.
#' @return One-row tibble.
#' @examples
#' recalc_r2_floor_from_correlations(r_y = c(-0.53, 0.61), r2 = 0.41)
#' @export
recalc_r2_floor_from_correlations <- function(r_y, r2 = NULL,
                                              r_y_digits = 2, r2_digits = 2) {
  k <- length(r_y)
  inputs <- setNames(lapply(r_y, interval_from_digits, r_y_digits),
                     paste0("r", seq_len(k)))
  fn <- function(...) max(unlist(list(...))^2)
  recomp <- propagate_intervals(fn, inputs)
  out <- recalc_result("C4: R^2 >= max_i r_Yi^2", r2,
                       reported_interval(r2, r2_digits), recomp)
  if (!is.null(r2)) {
    r2_upper <- r2 + 0.5 * 10^(-r2_digits)
    out$consistent <- r2_upper >= recomp[["lower"]]
  }
  out
}

#' Recalculate Delta F for a hierarchical R^2 increment
#'
#' Identity: \eqn{\Delta F = ((R^2_2 - R^2_1)/\Delta k) / ((1 - R^2_2)/(N - k_2 - 1))}.
#' Composes with \code{recalc_p_from_f()} to recalculate the p-value for the
#' increment on (Delta k, N - k_2 - 1) df.
#'
#' @param r2_full,r2_reduced Numeric.
#' @param delta_k Numeric. Number of predictors added.
#' @param n,k_full Numeric.
#' @param delta_f Numeric. Reported Delta F (optional).
#' @param r2_full_digits,r2_reduced_digits,delta_k_digits,n_digits,k_full_digits,delta_f_digits Integer.
#' @return One-row tibble.
#' @examples
#' recalc_delta_f(r2_full = 0.41, r2_reduced = 0.18, delta_k = 1,
#'               n = 1923, k_full = 2, delta_f = 750.1)
#' @export
recalc_delta_f <- function(r2_full, r2_reduced, delta_k, n, k_full,
                           delta_f = NULL,
                           r2_full_digits = 2, r2_reduced_digits = 2,
                           delta_k_digits = 0, n_digits = 0,
                           k_full_digits = 0, delta_f_digits = 2) {
  recomp <- propagate_intervals(
    fn = function(r2_full, r2_reduced, delta_k, n, k_full) {
      ((r2_full - r2_reduced) / delta_k) / ((1 - r2_full) / (n - k_full - 1))
    },
    inputs = list(r2_full = interval_from_digits(r2_full, r2_full_digits),
                  r2_reduced = interval_from_digits(r2_reduced,
                                                    r2_reduced_digits),
                  delta_k = interval_from_digits(delta_k, delta_k_digits),
                  n = interval_from_digits(n, n_digits),
                  k_full = interval_from_digits(k_full, k_full_digits))
  )
  recalc_result("C5: Delta F = ((R^2_2 - R^2_1)/Dk) / ((1 - R^2_2)/(N-k_2-1))",
                delta_f, reported_interval(delta_f, delta_f_digits), recomp)
}

#' Cross-check the beta-r decomposition against F (no R^2 needed)
#'
#' Combines B1 and B4: \eqn{\sum \beta_i r_{Yi} = Fk/(Fk + N - k - 1)}.
#' When the paper omits R^2 but reports F along with betas and r_Yi values,
#' this is the most economical check available.
#'
#' @param betas Numeric vector. Standardized regression coefficients.
#' @param r_y Numeric vector. Zero-order correlations with Y.
#' @param f,n,k Numeric.
#' @param betas_digits,r_y_digits,f_digits,n_digits,k_digits Integer.
#' @return One-row tibble.
#' @examples
#' recalc_betar_from_f(betas = c(-0.39, 0.47), r_y = c(-0.53, 0.61),
#'                    f = 665.4, n = 1923, k = 2)
#' @export
recalc_betar_from_f <- function(betas, r_y, f, n, k,
                                betas_digits = 2, r_y_digits = 2,
                                f_digits = 2, n_digits = 0, k_digits = 0) {
  stopifnot(length(betas) == length(r_y))
  k_pred <- length(betas)
  rep_int <- propagate_intervals(
    fn = function(f, n, k) (f * k) / (f * k + (n - k - 1)),
    inputs = list(f = interval_from_digits(f, f_digits),
                  n = interval_from_digits(n, n_digits),
                  k = interval_from_digits(k, k_digits))
  )
  inputs <- c(
    setNames(lapply(betas, interval_from_digits, betas_digits),
             paste0("b", seq_len(k_pred))),
    setNames(lapply(r_y, interval_from_digits, r_y_digits),
             paste0("r", seq_len(k_pred)))
  )
  fn <- function(...) {
    args <- list(...)
    bs <- unlist(args[paste0("b", seq_len(k_pred))])
    rs <- unlist(args[paste0("r", seq_len(k_pred))])
    sum(bs * rs)
  }
  recomp <- propagate_intervals(fn, inputs)
  recalc_result("C6: sum beta_i r_Yi = Fk/(Fk + N-k-1)",
                NA_real_,
                c(rep_int[["lower"]], rep_int[["upper"]]),
                recomp)
}
