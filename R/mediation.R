#' Sobel-family p-value for an indirect effect
#'
#' Recovers each path's SE by inverting the reported per-path p-value
#' (\eqn{SE_a = |a| / q(1 - p_a/2)} for two-tailed tests, where \eqn{q} is
#' the inverse z or t CDF), forms \eqn{ab} and the delta-method SE, and
#' returns the two-tailed p for the indirect effect.
#'
#' Three variants of \eqn{Var(ab)} are supported:
#' \describe{
#'   \item{\code{"sobel"}}{First-order delta:
#'     \eqn{Var(ab) = b^2 SE_a^2 + a^2 SE_b^2}.}
#'   \item{\code{"aroian"}}{Adds the second-order term:
#'     \eqn{+\, SE_a^2 SE_b^2}.}
#'   \item{\code{"goodman"}}{Subtracts it; clamped at 0 to avoid negative
#'     variance: \eqn{-\, SE_a^2 SE_b^2}.}
#' }
#'
#' Non-finite intermediates (e.g. when a reported \eqn{p_a} of exactly 1
#' makes \eqn{SE_a} infinite, or when both paths are reported as exactly 0)
#' return \code{NA_real_}.
#'
#' @param a,b Numeric. Path coefficients.
#' @param p_a,p_b Numeric. Reported p-values for paths a and b.
#' @param df Numeric scalar or \code{NULL}. If \code{NULL}, the reference
#'   distribution is standard normal; otherwise Student t with this df.
#' @param two_tailed Logical. Whether the reported per-path p-values are
#'   two-tailed (default \code{TRUE}).
#' @param method One of \code{"sobel"}, \code{"aroian"}, \code{"goodman"}.
#' @return Numeric scalar, the two-tailed p for the indirect effect.
#' @keywords internal
sobel_p <- function(
  a,
  b,
  p_a,
  p_b,
  df = NULL,
  two_tailed = TRUE,
  method = c("sobel", "aroian", "goodman")
) {
  method <- match.arg(method)
  tail_factor <- if (two_tailed) 2 else 1
  crit <- if (is.null(df)) {
    function(p) stats::qnorm(p)
  } else {
    function(p) stats::qt(p, df = df)
  }
  q_a <- crit(1 - p_a / tail_factor)
  q_b <- crit(1 - p_b / tail_factor)
  SE_a <- abs(a) / q_a
  SE_b <- abs(b) / q_b
  var_ab <- b^2 * SE_a^2 + a^2 * SE_b^2
  var_ab <- switch(
    method,
    sobel = var_ab,
    aroian = var_ab + SE_a^2 * SE_b^2,
    goodman = max(var_ab - SE_a^2 * SE_b^2, 0)
  )
  SE_ab <- sqrt(var_ab)
  ab <- a * b
  z <- ab / SE_ab
  if (!is.finite(z)) {
    return(NA_real_)
  }
  if (is.null(df)) {
    2 * stats::pnorm(-abs(z))
  } else {
    2 * stats::pt(-abs(z), df = df)
  }
}

# Input interval for a reported p-value: rounding interval if reported as
# equality, one-sided interval if reported as p < x or p > x. Clamped to
# (.Machine$double.eps, 1) so that inverse-CDF calls remain finite.
p_input_interval <- function(p, digits, op = "eq", rounding = "either") {
  eps <- .Machine$double.eps
  switch(
    op,
    eq = interval_from_digits(p, digits, lo = eps, hi = 1, rounding = rounding),
    lt = ,
    le = c(eps, p),
    gt = ,
    ge = c(p, 1),
    stop("Unknown op: ", op)
  )
}

# Require digit arguments to be supplied explicitly. NULL is the documented
# default and forces the caller to state the reporting precision rather than
# silently inheriting a guess.
require_digits <- function(...) {
  args <- list(...)
  missing <- vapply(args, is.null, logical(1))
  if (any(missing)) {
    stop(
      "Missing required digits argument(s): ",
      paste(names(args)[missing], collapse = ", "),
      ". Specify the number of decimal places each value was reported to.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Verify a reported indirect effect against the product of path coefficients
#'
#' Identity: \eqn{ab = a \cdot b}. Catches transcription errors in the
#' reported indirect effect by checking whether the rounding interval of
#' the reported \eqn{ab} overlaps the interval implied by the rounding
#' intervals of \eqn{a} and \eqn{b}.
#'
#' @param a,b Numeric. Path coefficients.
#' @param ab Numeric. Reported indirect effect (optional).
#' @param a_digits,b_digits,ab_digits Integer. Number of decimal places each
#'   value was reported to. No defaults: the caller must specify the
#'   precision of every value they supply. \code{ab_digits} is required only
#'   when \code{ab} is provided.
#' @return One-row tibble.
#' @examples
#' recalc_mediation_ab(a = 0.268, b = -0.279, ab = -0.075,
#'                     a_digits = 3, b_digits = 3, ab_digits = 3)
#' @export
#' @param rounding See \code{\link{recalc_rounding}} for the accepted values.
recalc_mediation_ab <- function(
  a,
  b,
  ab = NULL,
  a_digits = NULL,
  b_digits = NULL,
  ab_digits = NULL,
  rounding = "either"
) {
  require_digits(a_digits = a_digits, b_digits = b_digits)
  if (!is.null(ab)) {
    require_digits(ab_digits = ab_digits)
  }
  recomp <- propagate_intervals(
    fn = function(a, b) a * b,
    inputs = list(
      a = interval_from_digits(a, a_digits, rounding = rounding),
      b = interval_from_digits(b, b_digits, rounding = rounding)
    )
  )
  recalc_result(
    "M1: ab = a * b",
    ab,
    reported_interval(ab, ab_digits, rounding = rounding),
    recomp
  )
}

#' Recalculate the Sobel-family p-value for an indirect effect
#'
#' Identity: \eqn{p = 2(1 - F(|ab| / SE_{ab}))} with
#' \eqn{SE_{ab}^2 = b^2 SE_a^2 + a^2 SE_b^2} (Sobel; \code{method}
#' chooses Aroian or Goodman variants). Each path's SE is recovered by
#' inverting its reported p-value, so this is a consistency check on what
#' the original paper would have computed under Sobel - not an endorsement
#' over bootstrap or Monte-Carlo CIs, which are now the field standard for
#' indirect-effect inference.
#'
#' Bounds on the recalculated p-value are obtained by corner enumeration of
#' the joint rounding intervals of \eqn{a}, \eqn{p_a}, \eqn{b}, \eqn{p_b}
#' (via \code{propagate_intervals()}). For each requested reference
#' distribution a separate row is returned; the z-based row is always
#' included.
#'
#' @param a,b Numeric. Path coefficients.
#' @param p_a,p_b Numeric. Reported per-path p-values.
#' @param p Numeric. Reported indirect-effect p-value (optional).
#' @param a_digits,p_a_digits,b_digits,p_b_digits,p_digits Integer. Number
#'   of decimal places each value was reported to. No defaults: the caller
#'   must specify the precision of every value they supply. \code{p_digits}
#'   is required only when \code{p} is provided.
#' @param p_a_op,p_b_op,p_op One of \code{"eq"}, \code{"lt"}, \code{"le"},
#'   \code{"gt"}, \code{"ge"}. Use \code{"lt"} for reports of the form
#'   \code{p < .001}.
#' @param two_tailed Logical. Whether the reported per-path p-values are
#'   two-tailed (default \code{TRUE}).
#' @param df_t Numeric vector or \code{NULL}. If \code{NULL}, only the
#'   z-based row is returned. Otherwise the z row plus one t-based row per
#'   element of \code{df_t}.
#' @param method Character vector of methods to run, any subset of
#'   \code{c("sobel", "aroian", "goodman")}. \code{NULL} (the default) runs
#'   all three.
#' @param label Optional string included as a \code{label} column.
#' @return Tibble, one row per method \eqn{\times} reference distribution.
#' @examples
#' # Rutherford et al. 2017, HAM-D outcome - all three Sobel variants
#' recalc_mediation_p(a = 0.268, p_a = 0.038,
#'                    b = -0.279, p_b = 0.005,
#'                    p = 0.046, df_t = c(32, 47),
#'                    a_digits = 3, p_a_digits = 3,
#'                    b_digits = 3, p_b_digits = 3, p_digits = 3,
#'                    label = "HAM-D")
#' @export
#' @param rounding See \code{\link{recalc_rounding}} for the accepted values.
recalc_mediation_p <- function(
  a,
  p_a,
  b,
  p_b,
  p = NULL,
  a_digits = NULL,
  p_a_digits = NULL,
  b_digits = NULL,
  p_b_digits = NULL,
  p_digits = NULL,
  p_a_op = "eq",
  p_b_op = "eq",
  p_op = "eq",
  two_tailed = TRUE,
  df_t = NULL,
  method = NULL,
  label = NULL,
  rounding = "either"
) {
  require_digits(
    a_digits = a_digits,
    p_a_digits = p_a_digits,
    b_digits = b_digits,
    p_b_digits = p_b_digits
  )
  if (!is.null(p)) {
    require_digits(p_digits = p_digits)
  }
  all_methods <- c("sobel", "aroian", "goodman")
  methods <- if (is.null(method)) {
    all_methods
  } else {
    match.arg(method, all_methods, several.ok = TRUE)
  }
  inputs <- list(
    a = interval_from_digits(a, a_digits, rounding = rounding),
    b = interval_from_digits(b, b_digits, rounding = rounding),
    p_a = p_input_interval(p_a, p_a_digits, op = p_a_op, rounding = rounding),
    p_b = p_input_interval(p_b, p_b_digits, op = p_b_op, rounding = rounding)
  )
  rep_int <- reported_interval(p, p_digits, op = p_op, rounding = rounding)

  one_row <- function(df_val, m) {
    fn <- function(a, b, p_a, p_b) {
      sobel_p(a, b, p_a, p_b, df = df_val, two_tailed = two_tailed, method = m)
    }
    recomp <- propagate_intervals(fn, inputs)
    out <- recalc_result(
      sprintf("M2: Sobel p (%s)", m),
      p,
      rep_int,
      recomp
    )
    out$method <- m
    out$reference <- if (is.null(df_val)) "z" else sprintf("t(%g)", df_val)
    out$df_t <- if (is.null(df_val)) NA_real_ else as.numeric(df_val)
    out
  }

  rows <- list()
  for (m in methods) {
    rows[[length(rows) + 1L]] <- one_row(NULL, m)
    for (d in df_t) {
      rows[[length(rows) + 1L]] <- one_row(d, m)
    }
  }
  out <- dplyr::bind_rows(rows)
  if (!is.null(label)) {
    out$label <- label
  }
  out
}

#' Recalculate a mediation indirect effect and its Sobel-family p-value
#'
#' Convenience wrapper that binds the row from \code{recalc_mediation_ab()}
#' (the \eqn{ab = a \cdot b} identity check) with the rows from
#' \code{recalc_mediation_p()} (the Sobel p-value check across reference
#' distributions and variants). See those functions for the underlying
#' identities and caveats.
#'
#' @param a,b Numeric. Path coefficients.
#' @param p_a,p_b Numeric. Reported per-path p-values.
#' @param ab Numeric. Reported indirect effect (optional).
#' @param p Numeric. Reported indirect-effect p-value (optional).
#' @param a_digits,p_a_digits,b_digits,p_b_digits,ab_digits,p_digits Integer.
#'   Number of decimal places each value was reported to. No defaults: the
#'   caller must specify the precision of every value they supply.
#'   \code{ab_digits} is required only when \code{ab} is provided;
#'   \code{p_digits} only when \code{p} is provided.
#' @param p_a_op,p_b_op,p_op Comparison operator for the corresponding
#'   reported p-value; see \code{recalc_mediation_p()}.
#' @param two_tailed Logical. Whether the reported per-path p-values are
#'   two-tailed (default \code{TRUE}).
#' @param df_t Numeric vector or \code{NULL}. Reference distributions for
#'   the Sobel rows; see \code{recalc_mediation_p()}.
#' @param method Character vector of methods to run, any subset of
#'   \code{c("sobel", "aroian", "goodman")}. \code{NULL} (the default) runs
#'   all three.
#' @param label Optional string included as a \code{label} column.
#' @return Tibble: one M1 row plus one M2 row per method \eqn{\times}
#'   reference distribution.
#' @examples
#' # Rutherford et al. 2017, HAM-D outcome - all three Sobel variants
#' recalc_mediation(a = 0.268, p_a = 0.038,
#'                  b = -0.279, p_b = 0.005,
#'                  ab = -0.075, p = 0.046,
#'                  df_t = c(32, 47),
#'                  a_digits = 3, p_a_digits = 3,
#'                  b_digits = 3, p_b_digits = 3,
#'                  ab_digits = 3, p_digits = 3,
#'                  label = "HAM-D")
#' @export
#' @param rounding See \code{\link{recalc_rounding}} for the accepted values.
recalc_mediation <- function(
  a,
  p_a,
  b,
  p_b,
  ab = NULL,
  p = NULL,
  a_digits = NULL,
  p_a_digits = NULL,
  b_digits = NULL,
  p_b_digits = NULL,
  ab_digits = NULL,
  p_digits = NULL,
  p_a_op = "eq",
  p_b_op = "eq",
  p_op = "eq",
  two_tailed = TRUE,
  df_t = NULL,
  method = NULL,
  label = NULL,
  rounding = "either"
) {
  row_ab <- recalc_mediation_ab(
    a,
    b,
    ab,
    a_digits = a_digits,
    b_digits = b_digits,
    ab_digits = ab_digits,
    rounding = rounding
  )
  row_ab$method <- NA_character_
  row_ab$reference <- NA_character_
  row_ab$df_t <- NA_real_
  rows_p <- recalc_mediation_p(
    a = a,
    p_a = p_a,
    b = b,
    p_b = p_b,
    p = p,
    a_digits = a_digits,
    p_a_digits = p_a_digits,
    b_digits = b_digits,
    p_b_digits = p_b_digits,
    p_digits = p_digits,
    p_a_op = p_a_op,
    p_b_op = p_b_op,
    p_op = p_op,
    two_tailed = two_tailed,
    df_t = df_t,
    method = method,
    label = NULL,
    rounding = rounding
  )
  out <- dplyr::bind_rows(row_ab, rows_p)
  if (!is.null(label)) {
    out$label <- label
  }
  out
}
