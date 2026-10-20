# === Odds-ratio multiverse engine (used by recalc_chisq_p / recalc_oddsratio_p)
#
# The odds-ratio track is an OPT-IN companion to the chi-square p-value
# multiverse. It reuses the same candidate-table grid (exact counts, or the
# input-rounding reconstruction from reported proportions) built by
# `.build_chisq_grid()` and the same gap-aware per-method union machinery
# (`.method_intervals` / `.apply_union_default` in infrastructure.R).
#
# A 2x2 table has rows = groups, columns = outcome levels:
#        level1  level2
#   g1     a       b
#   g2     c       d
# and the "native" odds ratio is (a * d) / (b * c): the odds of level 1 in
# group 1 relative to group 2.
#
# Design notes (see also the roxygen on `recalc_chisq_p`):
#   * Score / likelihood-ratio / exact-conditional tests of OR != 1 coincide
#     exactly with the `pearson`, `likelihood_ratio` and `fisher` methods
#     already in the chi-square p track, so the only genuinely new OR p-value
#     here is the Wald test on log(OR).
#   * Continuity correction adds a constant (`or_cc`, default 0.5; set 1 for the
#     less-common +1 variant) to every cell, applied only on zero cells
#     (`or_cc_apply = "on_zero"`, default), always, or never. It changes the
#     Wald estimate / CI / p; the conditional-MLE and the exact (Cornfield /
#     mid-p) CIs are computed on the raw integer table and ignore it.
#   * Orientation "both" additionally carries the reciprocal (1/OR) so that a
#     reported OR expressed against the other reference group is not mistaken
#     for a non-reproducible value.

# --- single-table primitives -------------------------------------------------

#' Finite value or NA (maps +/-Inf and NaN to NA_real_).
#' @keywords internal
.finite_or_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.finite(x), x, NA_real_)
}

#' Apply the continuity-correction policy to a 2x2 table.
#'
#' Returns a list with the (possibly corrected) table and a label describing
#' what was done.
#' @keywords internal
.or_apply_cc <- function(tab, cc, cc_apply) {
  has_zero <- any(tab == 0)
  do_cc <- switch(
    cc_apply,
    never = FALSE,
    always = TRUE,
    on_zero = has_zero,
    stop("Unknown cc_apply: ", cc_apply)
  )
  if (do_cc && cc > 0) {
    list(tab = tab + cc, label = paste0("cc:+", cc))
  } else {
    list(tab = tab, label = "cc:none")
  }
}

#' Unconditional-MLE odds ratio = (a d) / (b c).
#' @keywords internal
.or_uncond <- function(tab) {
  .finite_or_na((tab[1, 1] * tab[2, 2]) / (tab[1, 2] * tab[2, 1]))
}

#' Conditional-MLE odds ratio (noncentral hypergeometric), via fisher.test.
#'
#' Always computed on the raw integer table (fisher.test requires integers).
#' @keywords internal
.or_cond <- function(tab) {
  est <- tryCatch(
    suppressWarnings(unname(fisher.test(tab)$estimate)),
    error = function(e) NA_real_
  )
  .finite_or_na(est)
}

#' Wald two-sided p for log(OR) != 0.
#' @keywords internal
.or_wald_p <- function(tab) {
  a <- tab[1, 1]
  b <- tab[1, 2]
  cc <- tab[2, 1]
  d <- tab[2, 2]
  if (any(c(a, b, cc, d) <= 0)) {
    return(NA_real_)
  }
  logor <- log((a * d) / (b * cc))
  se <- sqrt(1 / a + 1 / b + 1 / cc + 1 / d)
  if (!is.finite(se) || se == 0) {
    return(NA_real_)
  }
  2 * pnorm(-abs(logor / se))
}

#' Wald log(OR) confidence interval.
#' @keywords internal
.or_ci_wald <- function(tab, alpha) {
  a <- tab[1, 1]
  b <- tab[1, 2]
  cc <- tab[2, 1]
  d <- tab[2, 2]
  if (any(c(a, b, cc, d) <= 0)) {
    return(c(NA_real_, NA_real_))
  }
  logor <- log((a * d) / (b * cc))
  se <- sqrt(1 / a + 1 / b + 1 / cc + 1 / d)
  z <- qnorm(1 - alpha / 2)
  c(.finite_or_na(exp(logor - z * se)), .finite_or_na(exp(logor + z * se)))
}

#' Cornfield / exact-conditional CI, via fisher.test on the raw table.
#' @keywords internal
.or_ci_fisher <- function(tab, alpha) {
  ci <- tryCatch(
    suppressWarnings(
      fisher.test(tab, conf.level = 1 - alpha)$conf.int
    ),
    error = function(e) c(NA_real_, NA_real_)
  )
  c(.finite_or_na(ci[1]), .finite_or_na(ci[2]))
}

#' Fisher noncentral hypergeometric probabilities over the support of `a`.
#'
#' For a 2x2 with row totals r1, r2 and level-1 column total k, under odds
#' ratio `psi`. Returns a named numeric vector (names are the support values).
#' @keywords internal
.nchg_probs <- function(r1, r2, k, psi) {
  lo <- max(0L, k - r2)
  hi <- min(k, r1)
  x <- lo:hi
  if (!is.finite(psi) || psi <= 0) {
    # Degenerate: psi -> 0 puts all mass at lo, psi -> Inf at hi.
    w <- if (psi <= 0) {
      as.numeric(x == lo)
    } else {
      as.numeric(x == hi)
    }
    return(stats::setNames(w / sum(w), x))
  }
  lw <- lchoose(r1, x) + lchoose(r2, k - x) + x * log(psi)
  w <- exp(lw - max(lw))
  stats::setNames(w / sum(w), x)
}

#' Mid-p exact (equal-tailed) odds-ratio CI.
#'
#' Inverts the mid-p conditional (noncentral hypergeometric) tail functions:
#' the lower limit solves \eqn{P(X \ge a) - 0.5 P(X = a) = \alpha/2} and the
#' upper limit solves \eqn{P(X \le a) - 0.5 P(X = a) = \alpha/2}. Matches the
#' mid-p OR interval reported by StatXact / OpenEpi / epitools.
#' @keywords internal
.or_ci_midp <- function(tab, alpha) {
  r1 <- sum(tab[1, ])
  r2 <- sum(tab[2, ])
  k <- sum(tab[, 1])
  a <- tab[1, 1]
  lo <- max(0L, k - r2)
  hi <- min(k, r1)
  if (r1 == 0 || r2 == 0 || k == 0 || (r1 + r2 - k) == 0) {
    return(c(NA_real_, NA_real_))
  }

  midp_upper_tail <- function(psi) {
    # P(X >= a) - 0.5 P(X = a)
    pr <- .nchg_probs(r1, r2, k, psi)
    xs <- as.integer(names(pr))
    sum(pr[xs >= a]) - 0.5 * pr[xs == a]
  }
  midp_lower_tail <- function(psi) {
    # P(X <= a) - 0.5 P(X = a)
    pr <- .nchg_probs(r1, r2, k, psi)
    xs <- as.integer(names(pr))
    sum(pr[xs <= a]) - 0.5 * pr[xs == a]
  }

  target <- alpha / 2
  # Solve on the log(psi) scale over a wide bracket.
  solve_limit <- function(fn, lower_limit) {
    g <- function(lp) fn(exp(lp)) - target
    lo_lp <- -50
    hi_lp <- 50
    f_lo <- g(lo_lp)
    f_hi <- g(hi_lp)
    if (!is.finite(f_lo) || !is.finite(f_hi) || f_lo * f_hi > 0) {
      return(if (lower_limit) 0 else Inf)
    }
    exp(stats::uniroot(g, lower = lo_lp, upper = hi_lp)$root)
  }

  lower <- if (a == lo) 0 else solve_limit(midp_upper_tail, TRUE)
  upper <- if (a == hi) Inf else solve_limit(midp_lower_tail, FALSE)
  c(.finite_or_na(lower), .finite_or_na(upper))
}

#' Profile-likelihood odds-ratio CI, via confint() on a logistic glm.
#'
#' This is the default CI returned by \code{confint()} on a binomial
#' \code{glm} in R and by SAS \code{PROC LOGISTIC}'s \code{CLODDS=PL}, so it is
#' the natural CI source for any OR reported out of a logistic model.
#' A table with a zero cell is quasi-completely separated: the logistic fit
#' diverges and \code{confint()} returns \code{NA} for one limit and an
#' arbitrary value for the other. Such tables are passed to
#' \code{.or_ci_profile_lr()}, which inverts the likelihood-ratio test
#' directly and gives the finite limit and an infinite (\code{NA}) one.
#' @keywords internal
.or_ci_profile <- function(tab, alpha) {
  if (any(tab == 0)) {
    return(.or_ci_profile_lr(tab, alpha))
  }
  succ <- c(tab[1, 1], tab[2, 1])
  fail <- c(tab[1, 2], tab[2, 2])
  grp <- factor(c("g1", "g2"), levels = c("g1", "g2"))
  ci <- tryCatch(
    suppressWarnings(suppressMessages({
      fit <- stats::glm(
        cbind(succ, fail) ~ grp,
        family = stats::binomial()
      )
      stats::confint(fit, parm = "grpg2", level = 1 - alpha)
    })),
    error = function(e) c(NA_real_, NA_real_)
  )
  # The grpg2 coefficient is log-odds(level1 | g2) - log-odds(level1 | g1)
  # = -log(OR_native), so exp(-rev(ci)) maps it back to the native OR scale.
  c(.finite_or_na(exp(-ci[2])), .finite_or_na(exp(-ci[1])))
}

#' Profile-likelihood odds-ratio CI by direct inversion of the LR test.
#'
#' With row 2's log-odds \eqn{\beta_0} and row 1's \eqn{\beta_0 + \theta}
#' (\eqn{\theta} = log of the native OR), each limit solves
#' \eqn{2(\ell_{max} - \max_{\beta_0} \ell(\beta_0, \theta)) =
#' \chi^2_{1, 1 - \alpha}}. Used for tables with a zero cell, where the MLE of
#' \eqn{\theta} is infinite: a zero in \code{b} or \code{c} makes it
#' \eqn{+\infty} (upper limit infinite), a zero in \code{a} or \code{d}
#' \eqn{-\infty} (lower limit 0). An infinite limit is returned as \code{NA},
#' and a zero lower limit as 0, as for the other CI methods.
#' @keywords internal
.or_ci_profile_lr <- function(tab, alpha) {
  a <- tab[1, 1]
  b <- tab[1, 2]
  cc <- tab[2, 1]
  d <- tab[2, 2]
  n1 <- a + b
  n2 <- cc + d
  if (n1 == 0 || n2 == 0 || a + cc == 0 || b + d == 0) {
    return(c(NA_real_, NA_real_))
  }
  log1pexp <- function(x) ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
  # binomial log-likelihood at the observed proportion, with 0 log 0 = 0
  ll_sat <- function(x, n) {
    if (x == 0 || x == n) 0 else x * log(x / n) + (n - x) * log(1 - x / n)
  }
  ll_max <- ll_sat(a, n1) + ll_sat(cc, n2)
  ll_profile <- function(theta) {
    ll <- function(b0) {
      a * (b0 + theta) - n1 * log1pexp(b0 + theta) + cc * b0 - n2 * log1pexp(b0)
    }
    stats::optimize(
      ll,
      c(-abs(theta) - 40, abs(theta) + 40),
      maximum = TRUE
    )$objective
  }
  q <- stats::qchisq(1 - alpha, df = 1)
  g <- function(theta) 2 * (ll_max - ll_profile(theta)) - q

  # theta is searched on [-50, 50] (ORs from e^-50 to e^50). The deviance is
  # below q at the MLE and above it far from it, so there is one root on each
  # side of a finite MLE, and one root in all on the finite side of an
  # infinite MLE. A finite MLE only arises when this is called directly (e.g.
  # to check it against confint()), since .or_ci_profile() sends only tables
  # with a zero cell here.
  theta_hat <- if (b == 0 || cc == 0) {
    Inf
  } else if (a == 0 || d == 0) {
    -Inf
  } else {
    log((a * d) / (b * cc))
  }
  root <- function(from, to) {
    f_from <- g(from)
    f_to <- g(to)
    if (!is.finite(f_from) || !is.finite(f_to) || f_from * f_to > 0) {
      return(NA_real_)
    }
    exp(stats::uniroot(g, c(from, to), tol = 1e-10)$root)
  }
  lower <- if (theta_hat == -Inf) 0 else root(-50, min(theta_hat, 50))
  upper <- if (theta_hat == Inf) NA_real_ else root(max(theta_hat, -50), 50)
  c(.finite_or_na(lower), .finite_or_na(upper))
}


# --- multiverse assembly ------------------------------------------------------

#' Build the odds-ratio multiverse results frame over the candidate-table grid.
#'
#' One row per (table, estimator, CI method, orientation, cc policy). Carries
#' the OR point estimate for that estimator, the CI bounds for that CI method,
#' and the Wald log(OR) p for the table; orientation "reciprocal" holds 1/OR
#' with swapped-and-inverted CI bounds (the Wald p is symmetric).
#' @keywords internal
.multiverse_or_results <- function(
  grid,
  or_estimators,
  or_ci_methods,
  or_orientation,
  or_alpha,
  or_cc,
  or_cc_apply
) {
  orientations <- if (identical(or_orientation, "both")) {
    c("as_is", "reciprocal")
  } else {
    "as_is"
  }

  rows <- vector("list", 0L)

  for (ti in seq_along(grid$tables)) {
    tab <- grid$tables[[ti]]
    lbl <- grid$labels[ti]

    cc_res <- .or_apply_cc(tab, or_cc, or_cc_apply)
    tab_cc <- cc_res$tab
    cc_label <- cc_res$label

    est_vals <- list(
      uncond_mle = .or_uncond(tab_cc),
      cond_mle = .or_cond(tab)
    )
    ci_vals <- list(
      wald_logor = .or_ci_wald(tab_cc, or_alpha),
      cornfield_exact = .or_ci_fisher(tab, or_alpha),
      midp_exact = .or_ci_midp(tab, or_alpha),
      profile_likelihood = .or_ci_profile(tab, or_alpha)
    )
    p_wald <- .or_wald_p(tab_cc)

    for (estimator in or_estimators) {
      for (ci_method in or_ci_methods) {
        for (orient in orientations) {
          or_native <- est_vals[[estimator]]
          ci_native <- ci_vals[[ci_method]]
          if (identical(orient, "reciprocal")) {
            or_v <- .finite_or_na(1 / or_native)
            ci_v <- c(
              .finite_or_na(1 / ci_native[2]),
              .finite_or_na(1 / ci_native[1])
            )
          } else {
            or_v <- or_native
            ci_v <- ci_native
          }
          rows[[length(rows) + 1L]] <- data.frame(
            source = "odds_ratio",
            estimator = estimator,
            ci_method = ci_method,
            orientation = orient,
            cc_label = cc_label,
            input_adj_counts = lbl,
            or_unrounded = or_v,
            ci_lower_unrounded = ci_v[1],
            ci_upper_unrounded = ci_v[2],
            p_unrounded = p_wald,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  if (length(rows) == 0L) {
    return(data.frame(
      source = character(),
      estimator = character(),
      ci_method = character(),
      orientation = character(),
      cc_label = character(),
      input_adj_counts = character(),
      or_unrounded = numeric(),
      ci_lower_unrounded = numeric(),
      ci_upper_unrounded = numeric(),
      p_unrounded = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  out <- do.call(rbind, rows)
  out[order(out$or_unrounded), , drop = FALSE]
}


# --- reproduced summary -------------------------------------------------------

#' One-row reproduced summary for the odds-ratio track (hull decisions).
#'
#' The gap-aware union decision for the OR point estimate is applied on top of
#' this by `.apply_union_default`; the CI-bound reconciliation is reported as a
#' convex hull of the achievable bounds.
#' @keywords internal
.reproduced_or_summary <- function(
  or_results,
  or_num,
  or_ci_lower_num,
  or_ci_upper_num,
  or_digits,
  rounding
) {
  fin <- function(x) x[is.finite(x)]

  est <- if (
    !is.null(or_results) &&
      nrow(or_results) > 0 &&
      "or_unrounded" %in% names(or_results)
  ) {
    fin(or_results$or_unrounded)
  } else {
    numeric(0)
  }
  lo <- if (!is.null(or_results) && nrow(or_results) > 0) {
    fin(or_results$ci_lower_unrounded)
  } else {
    numeric(0)
  }
  hi <- if (!is.null(or_results) && nrow(or_results) > 0) {
    fin(or_results$ci_upper_unrounded)
  } else {
    numeric(0)
  }

  min_or <- if (length(est)) min(est) else NA_real_
  max_or <- if (length(est)) max(est) else NA_real_
  min_ci_lower <- if (length(lo)) min(lo) else NA_real_
  max_ci_lower <- if (length(lo)) max(lo) else NA_real_
  min_ci_upper <- if (length(hi)) min(hi) else NA_real_
  max_ci_upper <- if (length(hi)) max(hi) else NA_real_

  out <- tibble::tibble(
    odds_ratio = or_num,
    min_or = min_or,
    max_or = max_or,
    min_or_rounded = .round_for_report(min_or, or_digits, rounding),
    max_or_rounded = .round_for_report(max_or, or_digits, rounding),
    or_ci_lower = or_ci_lower_num,
    or_ci_upper = or_ci_upper_num,
    min_or_ci_lower = min_ci_lower,
    max_or_ci_lower = max_ci_lower,
    min_or_ci_upper = min_ci_upper,
    max_or_ci_upper = max_ci_upper
  )

  out$or_inbounds <- if (is.na(or_num)) {
    NA
  } else {
    dplyr::between(or_num, out$min_or_rounded, out$max_or_rounded)
  }
  out$or_ci_lower_inbounds <- if (is.na(or_ci_lower_num)) {
    NA
  } else {
    dplyr::between(
      or_ci_lower_num,
      .round_for_report(min_ci_lower, or_digits, rounding),
      .round_for_report(max_ci_lower, or_digits, rounding)
    )
  }
  out$or_ci_upper_inbounds <- if (is.na(or_ci_upper_num)) {
    NA
  } else {
    dplyr::between(
      or_ci_upper_num,
      .round_for_report(min_ci_upper, or_digits, rounding),
      .round_for_report(max_ci_upper, or_digits, rounding)
    )
  }

  out
}


#' Run the full odds-ratio track and return the pieces that `recalc_chisq_p`
#' folds into its result list.
#'
#' Returns a list with `reproduced` (one row, OR-prefixed columns incl. the
#' gap-aware union decision and the OR-Wald p union), `or_results` (long), and
#' `or_method_intervals` (per-method OR-estimate intervals).
#' @keywords internal
.run_or_track <- function(
  grid,
  odds_ratio,
  or_ci_lower,
  or_ci_upper,
  or_digits,
  or_estimators,
  or_ci_methods,
  or_orientation,
  or_alpha,
  or_cc,
  or_cc_apply,
  p_num,
  p_operator,
  p_digits,
  rounding
) {
  or_results <- .multiverse_or_results(
    grid,
    or_estimators,
    or_ci_methods,
    or_orientation,
    or_alpha,
    or_cc,
    or_cc_apply
  )

  or_num <- if (!is.null(odds_ratio)) as.numeric(odds_ratio) else NA_real_
  or_ci_lower_num <- if (!is.null(or_ci_lower)) {
    as.numeric(or_ci_lower)
  } else {
    NA_real_
  }
  or_ci_upper_num <- if (!is.null(or_ci_upper)) {
    as.numeric(or_ci_upper)
  } else {
    NA_real_
  }

  reproduced <- .reproduced_or_summary(
    or_results,
    or_num,
    or_ci_lower_num,
    or_ci_upper_num,
    or_digits,
    rounding
  )

  # Gap-aware union for the OR point estimate. A "method" is a combination of
  # estimator, orientation, and continuity-correction policy; the CI-only
  # choices do not move the estimate.
  est_intervals <- .method_intervals(
    or_results,
    "or_unrounded",
    c("estimator", "orientation", "cc_label"),
    or_digits,
    rounding
  )
  applied <- .apply_union_default(
    reproduced,
    est_intervals,
    or_num,
    "equals",
    or_digits,
    inbounds_col = "or_inbounds",
    value_prefix = "or"
  )
  reproduced <- applied$reproduced

  # `.apply_union_default` adds un-prefixed union descriptors; prefix them so
  # they do not collide with the chi-square p track's identically-named columns
  # when the two reproduced rows are column-bound.
  union_cols <- c(
    "n_method_intervals",
    "covered_width",
    "hull_width",
    "covered_fraction"
  )
  for (nm in intersect(union_cols, names(reproduced))) {
    names(reproduced)[names(reproduced) == nm] <- paste0("or_", nm)
  }

  # Gap-aware union for the OR-specific (Wald log-OR) p-value, kept as its own
  # decision so the chi-square association p stays untouched.
  p_intervals <- .method_intervals(
    or_results,
    "p_unrounded",
    c("orientation", "cc_label"),
    p_digits,
    rounding
  )
  p_union <- .union_summary(p_intervals, p_num, p_operator, p_digits)
  reproduced$or_p <- p_num
  reproduced$or_p_inbounds <- p_union$inbounds_union
  reproduced$or_p_covered_fraction <- p_union$covered_fraction

  list(
    reproduced = reproduced,
    or_results = or_results,
    or_method_intervals = applied$method_intervals
  )
}
