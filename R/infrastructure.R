#' Rounding interval implied by a value reported to `digits` decimal places
#'
#' Optional `lo`/`hi` clamping is useful for naturally bounded quantities
#' (e.g. `p` values reported to 3 dp near 0 or 1).
#' @keywords internal
interval_from_digits <- function(value, digits, lo = -Inf, hi = Inf) {
  delta <- 0.5 * 10^(-digits)
  c(max(value - delta, lo), min(value + delta, hi))
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
#' Use `op = "lt"` for reports of the form `p < .001`.
#' @keywords internal
reported_interval <- function(value, digits, op = "eq") {
  if (is.null(value) || is.na(value)) return(c(NA_real_, NA_real_))
  delta <- 0.5 * 10^(-digits)
  switch(
    op,
    eq = c(value - delta, value + delta),
    lt = c(-Inf, value),
    le = c(-Inf, value),
    gt = c(value, Inf),
    ge = c(value, Inf),
    stop("Unknown op: ", op)
  )
}

#' Standard return shape: tibble with reported and recomputed intervals
#' @keywords internal
recalc_result <- function(check, reported_value, reported_int, recomputed_int) {
  if (anyNA(reported_int)) {
    consistent <- NA
  } else {
    consistent <- (recomputed_int[["upper"]] >= reported_int[1]) &
                  (recomputed_int[["lower"]] <= reported_int[2])
  }
  tibble::tibble(
    check = check,
    reported = if (is.null(reported_value)) NA_real_ else as.numeric(reported_value),
    reported_lower = reported_int[1],
    reported_upper = reported_int[2],
    recomputed_lower = recomputed_int[["lower"]],
    recomputed_upper = recomputed_int[["upper"]],
    consistent = consistent
  )
}
