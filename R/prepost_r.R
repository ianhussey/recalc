#' Implied pre-post (within-person) correlation from reported SDs
#'
#' For a pre/post (repeated-measures) design the change-score variance is tied
#' to the pre and post SDs and the within-person correlation \eqn{r} by
#' \deqn{Var(change) = SD_{pre}^2 + SD_{post}^2 - 2 r\, SD_{pre} SD_{post},}
#' so a reported change SD pins down \eqn{r}:
#' \deqn{r = (SD_{pre}^2 + SD_{post}^2 - SD_{change}^2) / (2\, SD_{pre} SD_{post}).}
#'
#' This is the analytic inverse of the quantity DeBruine's \pkg{within}
#' package/app obtains by simulation (it builds an exact-statistics sample with
#' \code{faux::rnorm_multi(..., empirical = TRUE)} and solves for \eqn{r} by
#' optimisation); the closed form here avoids optimiser tolerance and returns
#' values outside \eqn{[-1, 1]} so impossibility is visible.
#'
#' \strong{Possibility.} A real change SD must lie in
#' \eqn{[\,|SD_{pre} - SD_{post}|,\ SD_{pre} + SD_{post}\,]} (the \eqn{r = +1}
#' and \eqn{r = -1} limits). When the whole recalculated interval falls above
#' \eqn{+1} or below \eqn{-1}, the three SDs cannot coexist. Whether a given
#' \eqn{r} is *plausible* (e.g. negative, or near 1) is a separate, measure-
#' dependent judgement left to the caller.
#'
#' Input rounding/truncation is propagated to the recalculated interval by
#' corner enumeration over the per-input preimage intervals (see
#' \code{\link{recalc_rounding}}); exact while the identity is monotone in each
#' input across the (sub-ULP) interval, which holds over typical reporting
#' ranges.
#'
#' @param sd_pre,sd_post,sd_change Numeric. Reported pre, post, and change-score
#'   SDs.
#' @param sd_pre_digits,sd_post_digits,sd_change_digits Integer. Decimal places
#'   each SD was reported to. Required.
#' @param r Numeric. A reported pre-post correlation to check (optional).
#' @param r_digits Integer. Decimal places \code{r} was reported to; required
#'   only when \code{r} is supplied.
#' @param rounding One of \code{"either"}, \code{"half_up"}, \code{"bankers"},
#'   \code{"truncate"}. See \code{\link{recalc_rounding}}.
#' @return One-row tibble (see \code{\link{recalc_result}}): the recalculated
#'   \code{[lower, upper]} interval for the implied \eqn{r}, and—if \code{r} is
#'   supplied—whether it is consistent with that interval.
#' @examples
#' # Gauhar (2016) TSC-40, change SD 13.01 with the TABLE component SDs:
#' recalc_prepost_r(sd_pre = 6.77, sd_post = 7.45, sd_change = 13.01,
#'                  sd_pre_digits = 2, sd_post_digits = 2, sd_change_digits = 2)
#' # ... and with the (incompatible) TEXT component SDs:
#' recalc_prepost_r(sd_pre = 14.58, sd_post = 7.79, sd_change = 13.01,
#'                  sd_pre_digits = 2, sd_post_digits = 2, sd_change_digits = 2)
#' @export
recalc_prepost_r <- function(sd_pre, sd_post, sd_change,
                             sd_pre_digits = NULL, sd_post_digits = NULL,
                             sd_change_digits = NULL,
                             r = NULL, r_digits = NULL, rounding = "either") {
  rounding <- match.arg(rounding, .recalc_rounding_choices)
  require_digits(sd_pre_digits = sd_pre_digits, sd_post_digits = sd_post_digits,
                 sd_change_digits = sd_change_digits)
  if (!is.null(r)) require_digits(r_digits = r_digits)

  inputs <- list(
    sd_pre    = interval_from_digits(sd_pre,    sd_pre_digits,    lo = 0, rounding = rounding),
    sd_post   = interval_from_digits(sd_post,   sd_post_digits,   lo = 0, rounding = rounding),
    sd_change = interval_from_digits(sd_change, sd_change_digits, lo = 0, rounding = rounding)
  )
  fn <- function(sd_pre, sd_post, sd_change) {
    (sd_pre^2 + sd_post^2 - sd_change^2) / (2 * sd_pre * sd_post)
  }
  recomp <- propagate_intervals(fn, inputs)
  recalc_result(
    "implied pre-post r = (sd_pre^2 + sd_post^2 - sd_change^2) / (2 sd_pre sd_post)",
    r, reported_interval(r, r_digits, rounding = rounding), recomp)
}

#' Change-score SD expected for a given pre-post correlation (reverse direction)
#'
#' \eqn{SD_{change} = \sqrt{SD_{pre}^2 + SD_{post}^2 - 2 r\, SD_{pre} SD_{post}}}.
#' The complement of \code{\link{recalc_prepost_r}}: given component SDs and a
#' correlation, what change SD should appear? Useful for checking a reported
#' change SD against the change SD a normal \eqn{r} (say 0.5) would produce.
#'
#' @inheritParams recalc_prepost_r
#' @param sd_pre,sd_post Numeric. Reported pre-test and post-test SDs.
#' @param sd_pre_digits,sd_post_digits Integer. Decimal places \code{sd_pre} /
#'   \code{sd_post} were reported to. Required.
#' @param r Numeric. Pre-post correlation to assume.
#' @param r_digits Integer. Decimal places \code{r} was reported to. Required.
#' @param sd_change Numeric. A reported change SD to check (optional).
#' @param sd_change_digits Integer. Required only when \code{sd_change} is given.
#' @return One-row tibble (see \code{\link{recalc_result}}).
#' @examples
#' recalc_change_sd_from_r(sd_pre = 14.58, sd_post = 7.79, r = 0.46,
#'                         sd_pre_digits = 2, sd_post_digits = 2, r_digits = 2)
#' @export
recalc_change_sd_from_r <- function(sd_pre, sd_post, r,
                                    sd_pre_digits = NULL, sd_post_digits = NULL,
                                    r_digits = NULL, sd_change = NULL,
                                    sd_change_digits = NULL, rounding = "either") {
  rounding <- match.arg(rounding, .recalc_rounding_choices)
  require_digits(sd_pre_digits = sd_pre_digits, sd_post_digits = sd_post_digits,
                 r_digits = r_digits)
  if (!is.null(sd_change)) require_digits(sd_change_digits = sd_change_digits)

  inputs <- list(
    sd_pre  = interval_from_digits(sd_pre,  sd_pre_digits,  lo = 0,  rounding = rounding),
    sd_post = interval_from_digits(sd_post, sd_post_digits, lo = 0,  rounding = rounding),
    r       = interval_from_digits(r,       r_digits, lo = -1, hi = 1, rounding = rounding)
  )
  fn <- function(sd_pre, sd_post, r) {
    sqrt(max(sd_pre^2 + sd_post^2 - 2 * r * sd_pre * sd_post, 0))
  }
  recomp <- propagate_intervals(fn, inputs)
  recalc_result(
    "change SD = sqrt(sd_pre^2 + sd_post^2 - 2 r sd_pre sd_post)",
    sd_change, reported_interval(sd_change, sd_change_digits, rounding = rounding), recomp)
}

#' Implied common pre-post correlation from a 2x2 RM-ANOVA interaction F
#'
#' For two timepoints the Group x Time interaction \eqn{F} of a 2 (group) x 2
#' (time) repeated-measures ANOVA equals \eqn{t^2} of an independent-samples
#' test on the change scores, so \eqn{F} fixes the pooled change-score SD
#' \deqn{SD_{pool} = |\bar{d}_1 - \bar{d}_2| / (\sqrt{F}\,\sqrt{1/n_1 + 1/n_2}),}
#' where \eqn{\bar{d}_g} is the group-\eqn{g} mean change. Assuming a single
#' common within-person \eqn{r} across groups, the pooled change variance is
#' linear in \eqn{r} and is solved in closed form. As with
#' \code{\link{recalc_prepost_r}}, an \eqn{r} outside \eqn{[-1, 1]} signals an
#' \eqn{F} incompatible with the reported means and SDs.
#'
#' @param f Numeric. Reported Group x Time interaction F.
#' @param m1b,sd1b,m1p,sd1p,n1 Group 1 baseline/post means and SDs, and n.
#' @param m2b,sd2b,m2p,sd2p,n2 Group 2 baseline/post means and SDs, and n.
#' @param f_digits,m_digits,sd_digits Integer. Decimal places F, the means, and
#'   the SDs were reported to. Required. \code{n1}, \code{n2} are treated as
#'   exact.
#' @param r,r_digits Optional reported r to check, and its decimal places.
#' @param rounding See \code{\link{recalc_rounding}}.
#' @return One-row tibble (see \code{\link{recalc_result}}).
#' @examples
#' # Pu et al. (2026) HAMD-17: F(1,76) = 117.055 -> r ~ 0.885
#' recalc_prepost_r_from_f(
#'   f = 117.055,
#'   m1b = 24.71, sd1b = 4.137, m1p = 8.96,  sd1p = 5.237, n1 = 52,
#'   m2b = 24.81, sd2b = 3.774, m2p = 15.95, sd2p = 5.714, n2 = 26,
#'   f_digits = 3, m_digits = 2, sd_digits = 3)
#' @export
recalc_prepost_r_from_f <- function(f, m1b, sd1b, m1p, sd1p, n1,
                                       m2b, sd2b, m2p, sd2p, n2,
                                    f_digits = NULL, m_digits = NULL,
                                    sd_digits = NULL, r = NULL, r_digits = NULL,
                                    rounding = "either") {
  rounding <- match.arg(rounding, .recalc_rounding_choices)
  require_digits(f_digits = f_digits, m_digits = m_digits, sd_digits = sd_digits)
  if (!is.null(r)) require_digits(r_digits = r_digits)

  mi  <- function(v) interval_from_digits(v, m_digits,  rounding = rounding)
  sdi <- function(v) interval_from_digits(v, sd_digits, lo = 0, rounding = rounding)
  inputs <- list(
    f   = interval_from_digits(f, f_digits, lo = 0, rounding = rounding),
    m1b = mi(m1b), sd1b = sdi(sd1b), m1p = mi(m1p), sd1p = sdi(sd1p),
    m2b = mi(m2b), sd2b = sdi(sd2b), m2p = mi(m2p), sd2p = sdi(sd2p)
  )
  fn <- function(f, m1b, sd1b, m1p, sd1p, m2b, sd2b, m2p, sd2p) {
    diff_change <- (m1p - m1b) - (m2p - m2b)
    s_pooled    <- abs(diff_change) / (sqrt(f) * sqrt(1 / n1 + 1 / n2))
    Nden <- n1 + n2 - 2
    A <- ((n1 - 1) * (sd1b^2 + sd1p^2) + (n2 - 1) * (sd2b^2 + sd2p^2)) / Nden
    B <- ((n1 - 1) * (2 * sd1b * sd1p) + (n2 - 1) * (2 * sd2b * sd2p)) / Nden
    (A - s_pooled^2) / B
  }
  recomp <- propagate_intervals(fn, inputs)
  recalc_result(
    "implied common pre-post r from 2x2 (group x time) RM-ANOVA interaction F",
    r, reported_interval(r, r_digits, rounding = rounding), recomp)
}

#' Recalculate partial eta-squared from an F statistic
#'
#' \eqn{\eta_p^2 = (F\, df_{effect}) / (F\, df_{effect} + df_{error})}. For a
#' 1-df effect this also pins the total N implied by the test
#' (\eqn{df_{error} = N - \#groups}), a useful cross-check against the reported
#' sample size.
#'
#' @param f Numeric. Reported F.
#' @param df_effect,df_error Numeric. Effect and error degrees of freedom
#'   (treated as exact).
#' @param f_digits Integer. Decimal places F was reported to. Required.
#' @param eta,eta_digits Optional reported partial eta-squared to check, and its
#'   decimal places.
#' @param rounding See \code{\link{recalc_rounding}}.
#' @return One-row tibble (see \code{\link{recalc_result}}).
#' @examples
#' recalc_partial_eta_from_f(f = 117.055, df_effect = 1, df_error = 76,
#'                           f_digits = 3, eta = 0.606, eta_digits = 3)
#' @export
recalc_partial_eta_from_f <- function(f, df_effect, df_error, f_digits = NULL,
                                      eta = NULL, eta_digits = NULL,
                                      rounding = "either") {
  rounding <- match.arg(rounding, .recalc_rounding_choices)
  require_digits(f_digits = f_digits)
  if (!is.null(eta)) require_digits(eta_digits = eta_digits)

  inputs <- list(f = interval_from_digits(f, f_digits, lo = 0, rounding = rounding))
  fn <- function(f) (f * df_effect) / (f * df_effect + df_error)
  recomp <- propagate_intervals(fn, inputs)
  recalc_result(
    "partial eta^2 = (F df_effect) / (F df_effect + df_error)",
    eta, reported_interval(eta, eta_digits, rounding = rounding), recomp)
}
