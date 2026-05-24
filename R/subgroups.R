#' Recalculate overall N, mean, and SD from subgroup statistics
#'
#' Identities:
#' \deqn{N = \sum_g n_g}
#' \deqn{M = \sum_g n_g M_g / N}
#' \deqn{\text{SD}^2 = \big(\sum_g (n_g - 1) s_g^2 + \sum_g n_g (M_g - M)^2\big)
#'                   / (N - 1)}
#'
#' The SD identity is the total-variance decomposition (within + between SS).
#' The textbook "pooled SD" formula
#' \eqn{\sqrt{\sum_g (n_g - 1) s_g^2 / \sum_g (n_g - 1)}} equals total SD only
#' when subgroup means coincide; using it as a check would false-flag any
#' paper whose subgroup means differ.
#'
#' Mismatch is a forensic signal: subgroup statistics were not computed from
#' the same observations as the overall statistics, a subgroup is missing or
#' misreported, or rounding was applied differently across rows.
#'
#' @param subgroup_ns Integer vector of subgroup sample sizes.
#' @param subgroup_means Numeric vector of subgroup means.
#' @param subgroup_sds Numeric vector of subgroup SDs (optional; if `NULL`
#'   the SD row is omitted).
#' @param overall_n,overall_mean,overall_sd Reported overall values (optional).
#' @param n_digits,mean_digits,sd_digits Integer. Number of decimal places each
#'   value was reported to. No defaults: the caller must specify the precision
#'   of every value they supply. `sd_digits` is required only when
#'   `subgroup_sds` is provided.
#' @return Tibble with one row per checked quantity.
#' @examples
#' recalc_total_from_subgroups(
#'   subgroup_ns    = c(16, 16),
#'   subgroup_means = c(4.50, 4.90),
#'   subgroup_sds   = c(2.30, 2.50),
#'   overall_n = 32, overall_mean = 4.70, overall_sd = 2.42,
#'   n_digits = 0, mean_digits = 2, sd_digits = 2
#' )
#' @export
recalc_total_from_subgroups <- function(subgroup_ns, subgroup_means,
                                        subgroup_sds = NULL,
                                        overall_n = NULL, overall_mean = NULL,
                                        overall_sd = NULL,
                                        n_digits = NULL, mean_digits = NULL,
                                        sd_digits = NULL) {
  require_digits(n_digits = n_digits, mean_digits = mean_digits)
  if (!is.null(subgroup_sds)) require_digits(sd_digits = sd_digits)

  k <- length(subgroup_ns)
  stopifnot(length(subgroup_means) == k, k >= 2)
  if (!is.null(subgroup_sds)) stopifnot(length(subgroup_sds) == k)

  n_names <- paste0("n", seq_len(k))
  m_names <- paste0("m", seq_len(k))
  inputs_nm <- c(
    setNames(lapply(subgroup_ns,    interval_from_digits, n_digits),    n_names),
    setNames(lapply(subgroup_means, interval_from_digits, mean_digits), m_names)
  )

  recomp_N <- propagate_intervals(
    fn = function(...) sum(unlist(list(...))),
    inputs = inputs_nm[n_names]
  )
  recomp_M <- propagate_intervals(
    fn = function(...) {
      args <- list(...)
      ns <- unlist(args[n_names]); ms <- unlist(args[m_names])
      sum(ns * ms) / sum(ns)
    },
    inputs = inputs_nm
  )

  out <- dplyr::bind_rows(
    recalc_result("E1: N = sum n_g",
                  overall_n, reported_interval(overall_n, n_digits), recomp_N),
    recalc_result("E1: M = (sum n_g M_g) / N",
                  overall_mean, reported_interval(overall_mean, mean_digits),
                  recomp_M)
  )

  if (!is.null(subgroup_sds)) {
    s_names <- paste0("s", seq_len(k))
    inputs_all <- c(
      inputs_nm,
      setNames(lapply(subgroup_sds, interval_from_digits, sd_digits), s_names)
    )
    recomp_SD <- propagate_intervals(
      fn = function(...) {
        args <- list(...)
        ns <- unlist(args[n_names])
        ms <- unlist(args[m_names])
        ss <- unlist(args[s_names])
        N <- sum(ns)
        M <- sum(ns * ms) / N
        within_ss  <- sum((ns - 1) * ss^2)
        between_ss <- sum(ns * (ms - M)^2)
        sqrt((within_ss + between_ss) / (N - 1))
      },
      inputs = inputs_all
    )
    out <- dplyr::bind_rows(out, recalc_result(
      "E1: SD = sqrt((sum (n_g-1) s_g^2 + sum n_g (M_g - M)^2) / (N - 1))",
      overall_sd, reported_interval(overall_sd, sd_digits), recomp_SD
    ))
  }
  out
}

#' Derive the implied N, mean, and SD of an unreported subgroup
#'
#' Given the overall N, M, and (optionally) SD and the reported subgroups,
#' invert the aggregation identities used by [recalc_total_from_subgroups()]
#' to recover what an unreported \eqn{k}-th subgroup would have to be for the
#' reported totals to be consistent.
#'
#' \deqn{n_\star = N - \sum_g n_g}
#' \deqn{M_\star = (N M - \sum_g n_g M_g) / n_\star}
#' \deqn{s_\star^2 = \big((N - 1) \text{SD}^2 - \sum_g (n_g - 1) s_g^2
#'                     - n_\star (M_\star - M)^2 - \sum_g n_g (M_g - M)^2\big)
#'                  / (n_\star - 1)}
#'
#' The `consistent` column reports physical-feasibility checks on the implied
#' missing-subgroup statistics:
#'
#' \describe{
#'   \item{`n_miss`}{`TRUE` iff the recomputed interval reaches \eqn{\geq 1}
#'         (i.e. some rounding makes the missing group non-empty); `FALSE`
#'         if the reported subgroup Ns already saturate or exceed
#'         `overall_n`.}
#'   \item{`M_miss`}{When `scale_min` and `scale_max` are supplied, `TRUE`
#'         iff the recomputed interval intersects `[scale_min, scale_max]`.
#'         `NA` when scale endpoints are not supplied (no bound to check).}
#'   \item{`SD_miss`}{`FALSE` if the propagated interval contains `NaN` —
#'         i.e. the implied variance is negative at some rounding corner, so
#'         no real subgroup reconciles the reported values. When `scale_min`
#'         and `scale_max` are supplied, additionally `FALSE` if the implied
#'         SD exceeds the maximum Bhatia-Davis bound
#'         \eqn{\sqrt{n_\star/(n_\star - 1)\,(M_\star - a)(b - M_\star)}}
#'         attainable over the rounding box (with `a = scale_min`,
#'         `b = scale_max`). The bound is the largest sample SD physically
#'         attainable by any \eqn{n_\star}-sample on \eqn{[a, b]} with mean
#'         \eqn{M_\star} (Bhatia & Davis, 2000), Bessel-corrected for the
#'         sample variance.}
#' }
#'
#' @param reported_ns,reported_means Numeric vectors. Statistics for the
#'   \eqn{k - 1} reported subgroups.
#' @param reported_sds Numeric vector (optional).
#' @param overall_n,overall_mean,overall_sd Numeric. Reported overall values.
#'   `overall_sd` is required to derive the SD row.
#' @param scale_min,scale_max Numeric. Optional logical endpoints of the
#'   measurement scale (e.g. `1` and `7` for a 7-point Likert). When supplied,
#'   the `M_miss` row is checked against `[scale_min, scale_max]` and the
#'   `SD_miss` row is additionally checked against the Bhatia-Davis bound. Both
#'   must be supplied together.
#' @param n_digits,mean_digits,sd_digits Integer. Number of decimal places each
#'   value was reported to. No defaults: the caller must specify the precision
#'   of every value they supply. `sd_digits` is required only when
#'   `reported_sds` and `overall_sd` are provided.
#' @return Tibble with one row per derived quantity. The `consistent` column
#'   reports the feasibility verdicts described in Details; `reported`
#'   columns remain `NA` (there is no reported missing-group value).
#' @references
#'   Bhatia, R., & Davis, C. (2000). A better bound on the variance.
#'   *American Mathematical Monthly*, 107(4), 353-357.
#' @examples
#' # 7-point Likert: bound checks via Bhatia-Davis and scale endpoints
#' recalc_missing_subgroup(
#'   reported_ns    = 50,
#'   reported_means = 3.50,
#'   reported_sds   = 1.40,
#'   overall_n = 100, overall_mean = 4.50, overall_sd = 1.50,
#'   scale_min = 1, scale_max = 7,
#'   n_digits = 0, mean_digits = 2, sd_digits = 2
#' )
#' @export
recalc_missing_subgroup <- function(reported_ns, reported_means,
                                    reported_sds = NULL,
                                    overall_n, overall_mean,
                                    overall_sd = NULL,
                                    scale_min = NULL, scale_max = NULL,
                                    n_digits = NULL, mean_digits = NULL,
                                    sd_digits = NULL) {
  require_digits(n_digits = n_digits, mean_digits = mean_digits)
  if (!is.null(reported_sds) && !is.null(overall_sd))
    require_digits(sd_digits = sd_digits)

  has_scale <- !is.null(scale_min) || !is.null(scale_max)
  if (has_scale) {
    if (is.null(scale_min) || is.null(scale_max))
      stop("scale_min and scale_max must both be supplied (or both NULL).")
    if (!(scale_min < scale_max))
      stop("scale_min must be strictly less than scale_max.")
  }

  k <- length(reported_ns)
  stopifnot(length(reported_means) == k, k >= 1)
  if (!is.null(reported_sds)) stopifnot(length(reported_sds) == k)

  n_names <- paste0("n", seq_len(k))
  m_names <- paste0("m", seq_len(k))
  inputs_nm <- c(
    list(N = interval_from_digits(overall_n,    n_digits),
         M = interval_from_digits(overall_mean, mean_digits)),
    setNames(lapply(reported_ns,    interval_from_digits, n_digits),    n_names),
    setNames(lapply(reported_means, interval_from_digits, mean_digits), m_names)
  )

  recomp_n <- propagate_intervals(
    fn = function(...) {
      args <- list(...)
      args$N - sum(unlist(args[n_names]))
    },
    inputs = inputs_nm[c("N", n_names)]
  )
  recomp_m <- propagate_intervals(
    fn = function(...) {
      args <- list(...)
      N <- args$N; M <- args$M
      ns <- unlist(args[n_names]); ms <- unlist(args[m_names])
      (N * M - sum(ns * ms)) / (N - sum(ns))
    },
    inputs = inputs_nm
  )

  out <- dplyr::bind_rows(
    recalc_result("E2: n_miss = N - sum n_g",
                  NULL, c(NA_real_, NA_real_), recomp_n),
    recalc_result("E2: M_miss = (N M - sum n_g M_g) / n_miss",
                  NULL, c(NA_real_, NA_real_), recomp_m)
  )

  # n_miss feasibility: must reach >= 1 at some rounding
  out$consistent[1] <- isTRUE(recomp_n[["upper"]] >= 1)
  # M_miss feasibility: intersects [scale_min, scale_max] (if scale supplied)
  if (has_scale) {
    out$consistent[2] <- isTRUE(
      recomp_m[["upper"]] >= scale_min &&
      recomp_m[["lower"]] <= scale_max
    )
  } else {
    out$consistent[2] <- NA
  }

  if (!is.null(reported_sds) && !is.null(overall_sd)) {
    s_names <- paste0("s", seq_len(k))
    inputs_all <- c(
      inputs_nm,
      list(SD = interval_from_digits(overall_sd, sd_digits)),
      setNames(lapply(reported_sds, interval_from_digits, sd_digits), s_names)
    )
    recomp_sd <- propagate_intervals(
      fn = function(...) {
        args <- list(...)
        N <- args$N; M <- args$M; SD <- args$SD
        ns <- unlist(args[n_names])
        ms <- unlist(args[m_names])
        ss <- unlist(args[s_names])
        n_miss <- N - sum(ns)
        m_miss <- (N * M - sum(ns * ms)) / n_miss
        var_miss <- ((N - 1) * SD^2 - sum((ns - 1) * ss^2)
                     - n_miss * (m_miss - M)^2 - sum(ns * (ms - M)^2)) /
                    (n_miss - 1)
        sqrt(var_miss)
      },
      inputs = inputs_all
    )
    out <- dplyr::bind_rows(out, recalc_result(
      "E2: SD_miss = sqrt(((N-1) SD^2 - within_g - between_g) / (n_miss - 1))",
      NULL, c(NA_real_, NA_real_), recomp_sd
    ))

    # SD_miss feasibility:
    #   - NaN -> FALSE (variance negative at some corner)
    #   - scale supplied -> additionally check against max Bhatia-Davis bound
    sd_lo <- recomp_sd[["lower"]]
    sd_hi <- recomp_sd[["upper"]]
    if (is.nan(sd_lo) || is.nan(sd_hi)) {
      out$consistent[3] <- FALSE
    } else if (has_scale) {
      # Max Bhatia-Davis bound over the rounding box. Corners where the
      # implied missing group is itself infeasible (n_miss < 2 or m_miss
      # outside the scale) contribute 0 and so do not raise the max.
      bd_recomp <- propagate_intervals(
        fn = function(...) {
          args <- list(...)
          N <- args$N; M <- args$M
          ns <- unlist(args[n_names]); ms <- unlist(args[m_names])
          n_miss <- N - sum(ns)
          m_miss <- (N * M - sum(ns * ms)) / n_miss
          if (!is.finite(n_miss) || n_miss < 2) return(0)
          if (!is.finite(m_miss) ||
              m_miss < scale_min || m_miss > scale_max) return(0)
          sqrt(n_miss / (n_miss - 1) *
               (m_miss - scale_min) * (scale_max - m_miss))
        },
        inputs = inputs_nm
      )
      bd_max <- bd_recomp[["upper"]]
      out$consistent[3] <- isTRUE(sd_lo <= bd_max)
    } else {
      out$consistent[3] <- TRUE
    }
  }
  out
}

#' Test for insufficient variance among subgroup SDs (Bartlett, lower tail)
#'
#' **Status: exploratory.** This function is offered as a candidate forensic
#' signal, not as a calibrated hypothesis test. The χ²(k−1) approximation
#' underlying Bartlett's K² is documented for the upper tail; lower-tail
#' calibration in the regimes psychology papers actually use (k = 2–5,
#' n_g = 20–100, mild non-normality) has not been validated. Use the output
#' as a flag for closer inspection of a table, not as evidence on its own.
#' See "Limitations" below.
#'
#' Bartlett's homogeneity-of-variance test, interpreted as a forensic check:
#' \deqn{K^2 = \big(\ln(s^2_\text{pooled}) \sum_g (n_g - 1)
#'                 - \sum_g (n_g - 1) \ln(s_g^2)\big) \,/\, C}
#' \deqn{C = 1 + \frac{1}{3(k - 1)}\left(\sum_g \frac{1}{n_g - 1}
#'                                       - \frac{1}{\sum_g (n_g - 1)}\right)}
#' Under a common-population null, \eqn{K^2 \sim \chi^2_{k-1}} approximately.
#' Subgroup SDs that cluster *more tightly* than sampling theory predicts
#' produce a small \eqn{K^2} and therefore a small *lower-tail* p-value. This
#' is TIVA's logic (Test of Insufficient Variance) applied at the
#' descriptive-statistic level rather than the test-statistic level.
#'
#' Sample SD has its own sampling variability — `SE(s)` is of order
#' \eqn{\sigma / \sqrt{2(n - 1)}} for normal data — and a paper whose
#' subgroup SDs do not exhibit that variability is a candidate for
#' interpolated or copy-pasted values.
#'
#' @section Limitations:
#' \itemize{
#'   \item **Calibration not validated.** Bartlett's χ²(k − 1) approximation
#'     was derived for upper-tail variance-heterogeneity testing. Its
#'     accuracy in the *lower* tail at small `n_g`, modest `k`, or under
#'     non-normality has not been established. The nominal `p < 0.05`
#'     threshold may not correspond to a 5\% false-positive rate.
#'   \item **Power against fabrication is unknown.** There is no published
#'     evidence that fabricated or interpolated SDs produce the specific
#'     pattern this test is sensitive to. Treat small p as a flag, not as a
#'     verdict.
#'   \item **Substantive heteroscedasticity is rare in honest data.** A
#'     small lower-tail p can arise from genuinely homoscedastic sampling
#'     (stratified designs, capped scales, treatment that homogenizes a
#'     dependent variable). The signal is "this looks unusual", not "this is
#'     fabricated".
#'   \item **Not yet wired into audit wrappers.** Intentionally excluded
#'     from any \code{audit_*()} dispatcher until the points above are
#'     resolved by simulation.
#' }
#'
#' @param subgroup_ns Integer vector.
#' @param subgroup_sds Numeric vector (same length as `subgroup_ns`, at least
#'   2 entries).
#' @param n_digits,sd_digits Integer. Number of decimal places each value was
#'   reported to. No defaults: the caller must specify the precision of every
#'   value they supply.
#' @return One-row tibble. `recomputed_lower` / `recomputed_upper` bracket the
#'   lower-tail p-value over the rounding box. Small values (e.g. < 0.05) are
#'   the forensic flag — subject to the caveats in "Limitations".
#' @examples
#' recalc_sd_concentration(subgroup_ns  = c(40, 40, 40, 40),
#'                         subgroup_sds = c(2.40, 2.40, 2.40, 2.40),
#'                         n_digits = 0, sd_digits = 2)
#' @export
recalc_sd_concentration <- function(subgroup_ns, subgroup_sds,
                                    n_digits = NULL, sd_digits = NULL) {
  require_digits(n_digits = n_digits, sd_digits = sd_digits)

  k <- length(subgroup_ns)
  stopifnot(length(subgroup_sds) == k, k >= 2)

  n_names <- paste0("n", seq_len(k))
  s_names <- paste0("s", seq_len(k))
  inputs <- c(
    setNames(lapply(subgroup_ns,  interval_from_digits, n_digits),  n_names),
    setNames(lapply(subgroup_sds, interval_from_digits, sd_digits), s_names)
  )

  fn <- function(...) {
    args <- list(...)
    ns <- unlist(args[n_names])
    ss <- unlist(args[s_names])
    df_g <- ns - 1
    df_tot <- sum(df_g)
    s2_pool <- sum(df_g * ss^2) / df_tot
    M <- df_tot * log(s2_pool) - sum(df_g * log(ss^2))
    C <- 1 + (sum(1 / df_g) - 1 / df_tot) / (3 * (k - 1))
    K2 <- M / C
    stats::pchisq(K2, df = k - 1)
  }
  recomp <- propagate_intervals(fn, inputs)
  recalc_result("E3: Bartlett K^2 lower-tail p (SD concentration)",
                NULL, c(NA_real_, NA_real_), recomp)
}
