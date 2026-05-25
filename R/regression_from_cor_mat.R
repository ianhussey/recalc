#' Convert the triangle of a correlation matrix to a full symmetric matrix
#'
#' Accepts either an upper or a lower triangle (the function auto-detects by
#' checking which half is all-NA) and returns the full symmetric matrix with
#' ones on the diagonal. Useful for converting a triangle reported in an
#' article into an input for `recalc_regression_from_cor_mat()`.
#'
#' @param triangle A square matrix or data frame. The unused half should
#'   contain `NA`.
#' @return A full symmetric correlation matrix, with row and column names
#'   preserved from the input.
#' @examples
#' lower <- matrix(NA, 3, 3,
#'                 dimnames = list(c("y", "x1", "x2"),
#'                                 c("y", "x1", "x2")))
#' lower[lower.tri(lower, diag = TRUE)] <- c(1, 0.5, 0.3, 1, 0.2, 1)
#' triangle_to_cor_matrix(lower)
#' @export
triangle_to_cor_matrix <- function(triangle) {
  if (is.data.frame(triangle)) triangle <- as.matrix(triangle)
  if (!is.matrix(triangle)) {
    stop("Input must be a matrix or data frame")
  }
  if (nrow(triangle) != ncol(triangle)) {
    stop("Input must be square")
  }
  upper_only <- all(is.na(triangle[lower.tri(triangle)])) &&
                all(!is.na(triangle[upper.tri(triangle)]))
  if (upper_only) triangle <- t(triangle)
  full <- matrix(0, nrow = nrow(triangle), ncol = ncol(triangle))
  full[lower.tri(full)] <- triangle[lower.tri(triangle)]
  full <- full + t(full)
  diag(full) <- 1
  dimnames(full) <- dimnames(triangle)
  full
}

#' Recalculate a linear regression from a reported correlation matrix
#'
#' Fits an OLS regression `Y ~ X1 + X2 + ...` directly from a reported
#' correlation matrix and sample size. Returns recalculated intervals on the
#' standardized regression coefficients, their Wald CIs, t statistics, and
#' p values. When a vector of SDs is supplied, unstandardized coefficients
#' \eqn{b_i = \beta_i \cdot \text{SD}_Y / \text{SD}_{X_i}} and their CIs are
#' also returned.
#'
#' Identities used (OLS from a standardized covariance matrix):
#' \eqn{\beta = R_{xx}^{-1} r_{xy}};
#' \eqn{R^2 = \beta^\top r_{xy}};
#' \eqn{\text{SE}(\beta_i) = \sqrt{(1-R^2)/(N-k-1) \cdot [R_{xx}^{-1}]_{ii}}};
#' \eqn{t_i = \beta_i / \text{SE}(\beta_i)};
#' \eqn{\text{CI}_{1-\alpha} = \beta_i \pm t_{1-\alpha/2,\,N-k-1} \cdot \text{SE}(\beta_i)}.
#'
#' Rounding intervals are propagated by enumerating
#' \eqn{2^{p(p-1)/2 + 1 + p_{\text{sd}}}} corner combinations, where \eqn{p} is
#' the number of variables in the model and \eqn{p_{\text{sd}}} is the number
#' of supplied SDs. Practical up to about \eqn{p = 6}.
#'
#' @param formula A Wilkinson-notation regression formula, e.g. `y ~ x1 + x2`.
#'   Variables must appear as row and column names of `cor_mat` (and of
#'   `sd_vector` when supplied).
#' @param cor_mat Numeric. A symmetric correlation matrix with row and column
#'   names. May contain extra variables not used in `formula`.
#' @param n Integer. Sample size.
#' @param sd_vector Optional named numeric vector of variable SDs. When
#'   supplied, unstandardized regression coefficients and their CIs are also
#'   recalculated. Must include all variables in `formula`.
#' @param reported_betas Optional named numeric vector of reported standardized
#'   regression coefficients, one per predictor. When supplied, each value is
#'   compared against the recalculated `Std.Beta` interval and `consistent` is
#'   filled in for those rows.
#' @param reported_b Optional named numeric vector of reported unstandardized
#'   regression coefficients, one per predictor. Compared against the
#'   recalculated `b` interval. Only meaningful when `sd_vector` is supplied
#'   (otherwise the function does not return `b` rows).
#' @param cor_mat_digits,n_digits,sd_vector_digits,reported_betas_digits,reported_b_digits Integer.
#'   Number of decimal places each value was reported to. `cor_mat_digits` has
#'   no default — the caller must state the precision of the reported
#'   correlations. `sd_vector_digits`, `reported_betas_digits`, and
#'   `reported_b_digits` are required only when their corresponding inputs
#'   are provided.
#' @param level Numeric. CI level. Default 0.95.
#' @param rounding See [recalc_rounding].
#' @return A long-format tibble with one row per (predictor x quantity). The
#'   `check` column names the identity recomputed (e.g. `"Std.Beta"`,
#'   `"Std.Beta CI lower"`, `"p"`, `"b"`); `predictor` is the predictor name
#'   from `formula`. `reported`, `reported_lower`, `reported_upper`, and
#'   `consistent` are populated for `Std.Beta` rows when `reported_betas` is
#'   supplied, and for `b` rows when `reported_b` is supplied; `NA` otherwise.
#' @examples
#' # Haller et al. (2022) correlation matrix and N
#' nm <- c("prosocial", "wellbeing", "support",
#'         "stress", "flexibility", "affect")
#' R <- matrix(c(1.00,  0.32,  0.34, -0.09,  0.21,  0.30,
#'               0.32,  1.00,  0.45, -0.54,  0.54,  0.64,
#'               0.34,  0.45,  1.00, -0.26,  0.32,  0.34,
#'              -0.09, -0.54, -0.26,  1.00, -0.52, -0.49,
#'               0.21,  0.54,  0.32, -0.52,  1.00,  0.50,
#'               0.30,  0.64,  0.34, -0.49,  0.50,  1.00),
#'             nrow = 6, dimnames = list(nm, nm))
#' recalc_regression_from_cor_mat(
#'   formula = prosocial ~ stress + affect + flexibility + support,
#'   cor_mat = R, n = 9496,
#'   cor_mat_digits = 2
#' )
#' @export
recalc_regression_from_cor_mat <- function(formula, cor_mat, n,
                                           sd_vector = NULL,
                                           reported_betas = NULL,
                                           reported_b = NULL,
                                           cor_mat_digits = NULL,
                                           n_digits = 0,
                                           sd_vector_digits = NULL,
                                           reported_betas_digits = NULL,
                                           reported_b_digits = NULL,
                                           level = 0.95,
                                           rounding = "either") {
  require_digits(cor_mat_digits = cor_mat_digits)
  if (!is.null(sd_vector)) require_digits(sd_vector_digits = sd_vector_digits)
  if (!is.null(reported_betas))
    require_digits(reported_betas_digits = reported_betas_digits)
  if (!is.null(reported_b))
    require_digits(reported_b_digits = reported_b_digits)
  stopifnot(inherits(formula, "formula"))
  stopifnot(is.matrix(cor_mat), nrow(cor_mat) == ncol(cor_mat))
  if (is.null(rownames(cor_mat)) || is.null(colnames(cor_mat))) {
    stop("cor_mat must have row and column names matching variables in formula")
  }
  if (!identical(rownames(cor_mat), colnames(cor_mat))) {
    stop("cor_mat row and column names must match")
  }

  vars <- all.vars(formula)
  y    <- vars[1]
  xs   <- vars[-1]
  k    <- length(xs)
  if (k < 1L) stop("formula must include at least one predictor")
  missing_vars <- setdiff(vars, rownames(cor_mat))
  if (length(missing_vars) > 0L) {
    stop("Variables not in cor_mat: ", paste(missing_vars, collapse = ", "))
  }

  R0 <- cor_mat[vars, vars, drop = FALSE]
  idx <- which(upper.tri(R0), arr.ind = TRUE)
  n_corr <- nrow(idx)
  corr_names <- sprintf("r_%d_%d", idx[, 1], idx[, 2])

  corr_intervals <- setNames(
    lapply(seq_len(n_corr), function(i) {
      interval_from_digits(R0[idx[i, 1], idx[i, 2]], cor_mat_digits,
                           lo = -1, hi = 1, rounding = rounding)
    }),
    corr_names
  )
  n_interval <- interval_from_digits(n, n_digits, lo = k + 2,
                                     rounding = rounding)

  use_sd <- !is.null(sd_vector)
  if (use_sd) {
    if (is.null(names(sd_vector))) {
      stop("sd_vector must be named to match variables in formula")
    }
    missing_sd <- setdiff(vars, names(sd_vector))
    if (length(missing_sd) > 0L) {
      stop("sd_vector missing variables: ",
           paste(missing_sd, collapse = ", "))
    }
    sd_ordered <- sd_vector[vars]
    sd_intervals <- setNames(
      lapply(seq_along(sd_ordered), function(i) {
        interval_from_digits(sd_ordered[i], sd_vector_digits, lo = 0,
                             rounding = rounding)
      }),
      paste0("sd_", seq_along(sd_ordered))
    )
  } else {
    sd_intervals <- NULL
  }

  validate_reported <- function(vec, label) {
    if (is.null(names(vec))) {
      stop(label, " must be a named numeric vector matching predictor names")
    }
    missing_p <- setdiff(xs, names(vec))
    if (length(missing_p) > 0L) {
      stop(label, " missing predictors: ",
           paste(missing_p, collapse = ", "))
    }
    vec[xs]
  }
  if (!is.null(reported_betas))
    reported_betas <- validate_reported(reported_betas, "reported_betas")
  if (!is.null(reported_b)) {
    if (!use_sd) {
      stop("reported_b requires sd_vector (b rows are only produced when ",
           "SDs are supplied)")
    }
    reported_b <- validate_reported(reported_b, "reported_b")
  }

  all_inputs <- c(corr_intervals, list(n = n_interval), sd_intervals)
  grid <- do.call(expand.grid,
                  c(all_inputs, list(KEEP.OUT.ATTRS = FALSE,
                                     stringsAsFactors = FALSE)))
  grid_mat <- as.matrix(grid)

  n_quant <- if (use_sd) 7L else 4L  # Std.Beta, CI_lo, CI_hi, p (+ b, b_CI_lo, b_CI_hi)
  n_out   <- k * n_quant
  vals    <- matrix(NA_real_, nrow = nrow(grid_mat), ncol = n_out)

  alpha <- 1 - level

  for (g in seq_len(nrow(grid_mat))) {
    row <- grid_mat[g, ]
    R <- diag(length(vars))
    for (i in seq_len(n_corr)) {
      v <- row[[corr_names[i]]]
      R[idx[i, 1], idx[i, 2]] <- v
      R[idx[i, 2], idx[i, 1]] <- v
    }
    n_val <- row[["n"]]
    Rxx <- R[-1, -1, drop = FALSE]
    rxy <- R[-1, 1]
    Rxx_inv <- tryCatch(solve(Rxx), error = function(e) NULL)
    if (is.null(Rxx_inv)) next
    beta <- as.numeric(Rxx_inv %*% rxy)
    r2 <- sum(beta * rxy)
    df_resid <- n_val - k - 1
    if (df_resid <= 0 || r2 >= 1 || r2 < 0) next
    diag_inv <- diag(Rxx_inv)
    if (any(diag_inv <= 0)) next
    se_beta <- sqrt((1 - r2) / df_resid * diag_inv)
    t_crit <- stats::qt(1 - alpha / 2, df_resid)
    ci_lo <- beta - t_crit * se_beta
    ci_hi <- beta + t_crit * se_beta
    t_val <- beta / se_beta
    p_val <- 2 * stats::pt(abs(t_val), df_resid, lower.tail = FALSE)

    out <- c(beta, ci_lo, ci_hi, p_val)
    if (use_sd) {
      sd_y <- row[["sd_1"]]
      sd_x <- vapply(seq_len(k),
                     function(i) row[[paste0("sd_", i + 1L)]],
                     numeric(1))
      b      <- beta  * sd_y / sd_x
      b_ci_l <- ci_lo * sd_y / sd_x
      b_ci_u <- ci_hi * sd_y / sd_x
      out <- c(out, b, b_ci_l, b_ci_u)
    }
    vals[g, ] <- out
  }

  lower <- suppressWarnings(apply(vals, 2, min, na.rm = TRUE))
  upper <- suppressWarnings(apply(vals, 2, max, na.rm = TRUE))
  lower[!is.finite(lower)] <- NA_real_
  upper[!is.finite(upper)] <- NA_real_

  # Each predictor contributes n_quant rows. For predictor i and quantity q,
  # the column index in vals is (q - 1) * k + i.
  std_checks <- c("Std.Beta", "Std.Beta CI lower",
                  "Std.Beta CI upper", "p")
  b_checks   <- c("b", "b CI lower", "b CI upper")
  quant_names <- if (use_sd) c(std_checks, b_checks) else std_checks

  rows <- vector("list", n_out)
  r <- 1L
  for (q in seq_along(quant_names)) {
    qname <- quant_names[q]
    for (i in seq_len(k)) {
      col <- (q - 1L) * k + i
      reported_value <- NA_real_
      reported_int <- c(NA_real_, NA_real_)
      if (qname == "Std.Beta" && !is.null(reported_betas)) {
        reported_value <- reported_betas[[i]]
        reported_int <- reported_interval(reported_value,
                                          reported_betas_digits,
                                          rounding = rounding)
      } else if (qname == "b" && !is.null(reported_b)) {
        reported_value <- reported_b[[i]]
        reported_int <- reported_interval(reported_value,
                                          reported_b_digits,
                                          rounding = rounding)
      }
      rows[[r]] <- recalc_result(
        check = qname,
        reported_value = reported_value,
        reported_int = reported_int,
        recalculated_int = c(lower = lower[col], upper = upper[col])
      )
      rows[[r]]$predictor <- xs[i]
      r <- r + 1L
    }
  }
  out <- dplyr::bind_rows(rows)
  out <- out[, c("check", "predictor", setdiff(names(out),
                                               c("check", "predictor")))]
  tibble::as_tibble(out)
}
