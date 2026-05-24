#' Rounding modes accepted package-wide
#'
#' Functions in `{recalc}` that propagate input rounding accept a `rounding`
#' argument with one of the following values:
#'
#' * `"either"` (default) — assume the reported value came from either
#'   round-half-up *or* bankers (round-half-to-even) rounding. Both produce
#'   the same closed preimage interval `[value − ½ULP, value + ½ULP]`, since
#'   the two methods agree everywhere except on the measure-zero tie set.
#' * `"half_up"` — same closed interval. Use when you want to be explicit
#'   that the source used round-half-up specifically.
#' * `"bankers"` — same closed interval. Use when you want to be explicit
#'   that the source used round-half-to-even specifically.
#' * `"truncate"` — assume the reported value was truncated toward zero.
#'   The preimage interval is **asymmetric**: `[value, value + ULP]` for
#'   `value > 0`, `[value − ULP, value]` for `value < 0`, and
#'   `[−ULP, +ULP]` for `value == 0`.
#'
#' Where `ULP = 10^(-digits)` is the unit of least precision implied by
#' the digit count.
#' @name recalc_rounding
#' @keywords internal
NULL

# Internal: valid rounding modes.
.recalc_rounding_choices <- c("either", "half_up", "bankers", "truncate")

# Internal: compute lower/upper offsets of the preimage interval relative to
# the reported value, for the given rounding mode. Returns c(lower, upper).
.rounding_offsets <- function(value, ulp, rounding) {
  if (rounding == "truncate") {
    if (value > 0)       c(0, ulp)
    else if (value < 0)  c(-ulp, 0)
    else                 c(-ulp, ulp)
  } else {
    # half_up, bankers, either: closed symmetric interval ±½ ULP
    c(-0.5 * ulp, 0.5 * ulp)
  }
}

#' Rounding interval implied by a value reported to `digits` decimal places
#'
#' Optional `lo`/`hi` clamping is useful for naturally bounded quantities
#' (e.g. `p` values reported to 3 dp near 0 or 1).
#'
#' @param value Numeric. The reported value.
#' @param digits Integer. Number of decimal places it was reported to.
#' @param lo,hi Numeric. Clamping bounds for naturally bounded quantities.
#' @param rounding One of `"either"`, `"half_up"`, `"bankers"`, `"truncate"`.
#'   See [recalc_rounding].
#' @keywords internal
interval_from_digits <- function(value, digits, lo = -Inf, hi = Inf,
                                 rounding = "either") {
  rounding <- match.arg(rounding, .recalc_rounding_choices)
  ulp <- 10^(-digits)
  off <- .rounding_offsets(value, ulp, rounding)
  c(max(value + off[1], lo), min(value + off[2], hi))
}

#' Propagate input rounding intervals through a function by corner enumeration
#'
#' Returns the min and max of `fn` over the Cartesian product of input
#' intervals. Exact when `fn` is monotone in each input over the interval
#' (which holds for every identity in this package over typical reporting
#' ranges).
#' @keywords internal
propagate_intervals <- function(fn, inputs) {
  grid <- do.call(expand.grid, c(inputs, list(KEEP.OUT.ATTRS = FALSE,
                                              stringsAsFactors = FALSE)))
  vals <- vapply(seq_len(nrow(grid)),
                 function(i) do.call(fn, as.list(grid[i, , drop = FALSE])),
                 numeric(1))
  c(lower = min(vals), upper = max(vals))
}

#' Reported interval implied by a value, its digits, and the comparison op
#'
#' Use `op = "lt"` for reports of the form `p < .001`. The `rounding`
#' argument is honored only when `op == "eq"`; for the inequality operators,
#' the threshold is taken at face value.
#'
#' @param value Numeric. The reported value.
#' @param digits Integer. Number of decimal places it was reported to.
#' @param op One of `"eq"`, `"lt"`, `"le"`, `"gt"`, `"ge"`.
#' @param rounding One of `"either"`, `"half_up"`, `"bankers"`, `"truncate"`.
#'   See [recalc_rounding].
#' @keywords internal
reported_interval <- function(value, digits, op = "eq", rounding = "either") {
  if (is.null(value) || is.na(value)) return(c(NA_real_, NA_real_))
  rounding <- match.arg(rounding, .recalc_rounding_choices)
  switch(
    op,
    eq = {
      ulp <- 10^(-digits)
      off <- .rounding_offsets(value, ulp, rounding)
      c(value + off[1], value + off[2])
    },
    lt = c(-Inf, value),
    le = c(-Inf, value),
    gt = c(value, Inf),
    ge = c(value, Inf),
    stop("Unknown op: ", op)
  )
}

#' Standard return shape: tibble with reported and recalculated intervals
#' @keywords internal
recalc_result <- function(check, reported_value, reported_int, recalculated_int) {
  if (anyNA(reported_int)) {
    consistent <- NA
  } else {
    consistent <- (recalculated_int[["upper"]] >= reported_int[1]) &
                  (recalculated_int[["lower"]] <= reported_int[2])
  }
  tibble::tibble(
    check = check,
    reported = if (is.null(reported_value)) NA_real_ else as.numeric(reported_value),
    reported_lower = reported_int[1],
    reported_upper = reported_int[2],
    recalculated_lower = recalculated_int[["lower"]],
    recalculated_upper = recalculated_int[["upper"]],
    consistent = consistent
  )
}
