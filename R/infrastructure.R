#' Rounding modes accepted package-wide
#'
#' Functions in `{recalc}` that propagate input rounding accept a `rounding`
#' argument with one of the following values:
#'
#' * `"either"` (default) - assume the reported value came from either
#'   round-half-up *or* bankers (round-half-to-even) rounding. Both produce
#'   the same closed preimage interval `[value - 1/2ULP, value + 1/2ULP]`, since
#'   the two methods agree everywhere except on the measure-zero tie set.
#' * `"half_up"` - same closed interval. Use when you want to be explicit
#'   that the source used round-half-up specifically.
#' * `"bankers"` - same closed interval. Use when you want to be explicit
#'   that the source used round-half-to-even specifically.
#' * `"truncate"` - assume the reported value was truncated toward zero.
#'   The preimage interval is **asymmetric**: `[value, value + ULP]` for
#'   `value > 0`, `[value - ULP, value]` for `value < 0`, and
#'   `[-ULP, +ULP]` for `value == 0`.
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
    if (value > 0) {
      c(0, ulp)
    } else if (value < 0) {
      c(-ulp, 0)
    } else {
      c(-ulp, ulp)
    }
  } else {
    # half_up, bankers, either: closed symmetric interval +/-1/2 ULP
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
interval_from_digits <- function(
  value,
  digits,
  lo = -Inf,
  hi = Inf,
  rounding = "either"
) {
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
  grid <- do.call(
    expand.grid,
    c(inputs, list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
  )
  vals <- vapply(
    seq_len(nrow(grid)),
    function(i) do.call(fn, as.list(grid[i, , drop = FALSE])),
    numeric(1)
  )
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
  if (is.null(value) || is.na(value)) {
    return(c(NA_real_, NA_real_))
  }
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

# === Shared: per-method intervals and gap-aware union ========================
#
# The recalculated values in a multiverse are not a continuous band: each
# distinct combination of *discrete* analytic choices (the test method, the
# effect-size formula, etc.) produces its own value, and across the *continuous*
# input-rounding variation each becomes a narrow interval. Between those
# per-method intervals there can be gaps that the global [min, max] hull fills
# in. The union of the per-method intervals is a subset of that hull but still
# contains every legitimately-producible value, so a union-based reproducibility
# decision is strictly sharper than a hull-based one at no cost in coverage.
# These helpers are shared by the chi-square, t-test p, and Cohen's d functions.

#' Per-"method" intervals over a multiverse results frame.
#'
#' Groups `results` by the discrete analytic-choice columns `group_cols` and,
#' for each group, returns the `[min, max]` of `value_col` across the remaining
#' (input-rounding) variation, plus those bounds rounded to reporting precision.
#' Returns a data frame with columns `method`, `v_min`, `v_max`,
#' `v_min_rounded`, `v_max_rounded` (callers rename the `v_*` columns), or NULL
#' if there is nothing to summarise.
#' @keywords internal
.method_intervals <- function(
  results,
  value_col,
  group_cols,
  digits,
  rounding
) {
  if (
    is.null(results) || nrow(results) == 0 || !value_col %in% names(results)
  ) {
    return(NULL)
  }
  group_cols <- intersect(group_cols, names(results))
  v <- results[[value_col]]
  keep <- is.finite(v)
  if (!any(keep)) {
    return(NULL)
  }
  v <- v[keep]
  if (length(group_cols) == 0) {
    grp <- factor(rep("all", length(v)))
  } else {
    grp <- interaction(
      results[keep, group_cols, drop = FALSE],
      drop = TRUE,
      sep = "|"
    )
  }
  vmin <- tapply(v, grp, min)
  vmax <- tapply(v, grp, max)
  out <- data.frame(
    method = names(vmin),
    v_min = as.numeric(vmin),
    v_max = as.numeric(vmax),
    stringsAsFactors = FALSE
  )
  out$v_min_rounded <- .round_for_report(out$v_min, digits, rounding)
  out$v_max_rounded <- .round_for_report(out$v_max, digits, rounding)
  out[order(out$v_min, out$v_max), , drop = FALSE]
}

#' Gap-aware union summary from per-method intervals.
#'
#' Merges the per-method rounded intervals (treating a gap of at most one
#' reporting ULP as contiguous) into disjoint intervals, then reports the
#' union-based `inbounds` decision plus the total covered width, the hull
#' width, and `covered_fraction = covered_width / hull_width` (small means the
#' hull is mostly gap). Widths are measured at reporting resolution: each
#' disjoint interval also covers its own ULP cell. For the inequality operators
#' the union adds no sharpness (a half-line), so the hull rule is used.
#' @keywords internal
.union_summary <- function(method_intervals, value_num, operator, digits) {
  ulp <- 10^(-digits)
  tol <- ulp * 1e-6
  if (is.null(method_intervals) || nrow(method_intervals) == 0) {
    return(tibble::tibble(
      inbounds_union = NA,
      n_method_intervals = NA_integer_,
      covered_width = NA_real_,
      hull_width = NA_real_,
      covered_fraction = NA_real_
    ))
  }
  lo <- method_intervals$v_min_rounded
  hi <- method_intervals$v_max_rounded
  ord <- order(lo, hi)
  lo <- lo[ord]
  hi <- hi[ord]

  m_lo <- lo[1]
  m_hi <- hi[1]
  merged_lo <- numeric(0)
  merged_hi <- numeric(0)
  if (length(lo) > 1) {
    for (i in 2:length(lo)) {
      if (lo[i] <= m_hi + ulp * 1.5) {
        m_hi <- max(m_hi, hi[i])
      } else {
        merged_lo <- c(merged_lo, m_lo)
        merged_hi <- c(merged_hi, m_hi)
        m_lo <- lo[i]
        m_hi <- hi[i]
      }
    }
  }
  merged_lo <- c(merged_lo, m_lo)
  merged_hi <- c(merged_hi, m_hi)

  covered_width <- sum((merged_hi - merged_lo) + ulp)
  hull_width <- (max(hi) - min(lo)) + ulp
  covered_fraction <- if (hull_width > 0) {
    covered_width / hull_width
  } else {
    NA_real_
  }

  if (is.na(value_num)) {
    inb <- NA
  } else if (operator == "equals") {
    inb <- any(value_num >= merged_lo - tol & value_num <= merged_hi + tol)
  } else if (operator %in% c("less_than", "less_than_or_equal_to")) {
    inb <- (min(lo) - tol) <= value_num
  } else {
    # greater_than(_or_equal_to)
    inb <- (max(hi) + tol) >= value_num
  }

  tibble::tibble(
    inbounds_union = inb,
    n_method_intervals = length(merged_lo),
    covered_width = covered_width,
    hull_width = hull_width,
    covered_fraction = covered_fraction
  )
}

#' Attach the union (default) decision + hull fallback + descriptors to a
#' one-row reproduced summary, and rename a method-intervals frame's value
#' columns. `inbounds_col` is e.g. "p_inbounds" or "d_inbounds".
#' @keywords internal
.apply_union_default <- function(
  reproduced,
  method_intervals,
  value_num,
  operator,
  digits,
  inbounds_col,
  value_prefix
) {
  union_out <- .union_summary(method_intervals, value_num, operator, digits)
  hull_col <- paste0(inbounds_col, "_hull")
  reproduced[[hull_col]] <- reproduced[[inbounds_col]] # keep hull for legibility
  reproduced[[inbounds_col]] <- union_out$inbounds_union # union is the default
  reproduced <- dplyr::bind_cols(
    reproduced,
    union_out[c(
      "n_method_intervals",
      "covered_width",
      "hull_width",
      "covered_fraction"
    )]
  )
  if (!is.null(method_intervals)) {
    names(method_intervals) <- sub(
      "^v_",
      paste0(value_prefix, "_"),
      names(method_intervals)
    )
  }
  list(reproduced = reproduced, method_intervals = method_intervals)
}

#' Standard return shape: tibble with reported and recalculated intervals
#' @keywords internal
recalc_result <- function(
  check,
  reported_value,
  reported_int,
  recalculated_int
) {
  if (anyNA(reported_int)) {
    consistent <- NA
  } else {
    consistent <- (recalculated_int[["upper"]] >= reported_int[1]) &
      (recalculated_int[["lower"]] <= reported_int[2])
  }
  tibble::tibble(
    check = check,
    reported = if (is.null(reported_value)) {
      NA_real_
    } else {
      as.numeric(reported_value)
    },
    reported_lower = reported_int[1],
    reported_upper = reported_int[2],
    recalculated_lower = recalculated_int[["lower"]],
    recalculated_upper = recalculated_int[["upper"]],
    consistent = consistent
  )
}
