#' Run every applicable recalc check on a regression report
#'
#' Takes the full set of quantities a paper might report and dispatches every
#' check whose inputs are all available. Per-predictor checks (A1-A4, A6 bound,
#' A8) run once per predictor whose inputs are present; model-level (B1-B5) and
#' cross-equation (C1, C2, C4, C6) checks run once when their inputs are
#' present. Missing inputs cause the corresponding checks to be skipped
#' silently. By default also runs the section-E label-confusion diagnostics
#' (E1 beta vs b, E2 R^2 vs adj R^2, E3 one- vs two-tailed p) where their
#' inputs are available.
#'
#' C5 (hierarchical Delta F) is structurally distinct (requires a nested-model
#' pair with separate R^2, k, and increment metadata) and is not dispatched by
#' this wrapper — call \code{recalc_delta_f()} directly.
#'
#' @param n,k Numeric scalars.
#' @param r2,adj_r2,f Numeric scalars (optional).
#' @param p_model Numeric scalar (optional). Reported model p.
#' @param p_model_op Comparison op for \code{p_model}.
#' @param sd_y Numeric scalar (optional).
#' @param b,se,t,p,beta Numeric vectors of length k (optional).
#' @param p_op Comparison op for per-predictor p.
#' @param ci_lower,ci_upper Numeric vectors of length k (optional).
#' @param r_y,sd_x,sr2 Numeric vectors of length k (optional).
#' @param predictor_names Character vector of length k (optional).
#' @param r_12 Numeric scalar (optional, k = 2 only).
#' @param delta_r2_blocks Numeric vector (optional, B5 only).
#' @param level Numeric. CI level. Default 0.95.
#' @param diagnostics Logical. If TRUE (default), run E1, E2, E3 diagnostics
#'   where their inputs are available.
#' @param n_digits,k_digits Integer.
#' @param r2_digits,adj_r2_digits,f_digits,p_model_digits Integer.
#' @param sd_y_digits Integer.
#' @param b_digits,se_digits,t_digits,p_digits Integer.
#' @param beta_digits,ci_digits Integer.
#' @param r_y_digits,sd_x_digits,sr2_digits Integer.
#' @param r_12_digits,delta_r2_blocks_digits Integer.
#' @return Tibble of every applicable check, with an extra \code{predictor}
#'   column (\code{"model"} or a predictor name/index).
#' @examples
#' audit_regression(
#'   n = 1923, k = 2, r2 = 0.41,
#'   beta = c(-0.39, 0.47), r_y = c(-0.53, 0.61),
#'   sd_y = 1.4, sd_x = c(1.2, 0.5),
#'   sd_y_digits = 1, sd_x_digits = 1,
#'   predictor_names = c("AI_reliance", "override_freq")
#' )
#' @export
audit_regression <- function(
    # Model-level
    n = NULL, k = NULL,
    r2 = NULL, adj_r2 = NULL, f = NULL,
    p_model = NULL, p_model_op = "eq",
    sd_y = NULL,
    # Per-predictor
    b = NULL, se = NULL, t = NULL, p = NULL, p_op = "eq",
    beta = NULL,
    ci_lower = NULL, ci_upper = NULL,
    r_y = NULL, sd_x = NULL, sr2 = NULL,
    predictor_names = NULL,
    # Cross-equation / hierarchical
    r_12 = NULL, delta_r2_blocks = NULL,
    level = 0.95,
    diagnostics = TRUE,
    # Digits
    n_digits = 0, k_digits = 0,
    r2_digits = 2, adj_r2_digits = 2, f_digits = 2, p_model_digits = 3,
    sd_y_digits = 2,
    b_digits = 2, se_digits = 2, t_digits = 2, p_digits = 3,
    beta_digits = 2, ci_digits = 2,
    r_y_digits = 2, sd_x_digits = 2, sr2_digits = 3,
    r_12_digits = 2, delta_r2_blocks_digits = 2,
    rounding = "either") {

  has <- function(x) !is.null(x) && length(x) >= 1L && !any(is.na(x))
  vget <- function(v, i) if (length(v) >= i) v[i] else NULL
  results <- list()
  add <- function(x, predictor = "model") {
    x$predictor <- predictor
    results[[length(results) + 1L]] <<- x
  }

  df <- if (has(n) && has(k)) n - k - 1 else NA_real_
  k_pred <- max(c(0L,
    length(b), length(se), length(t), length(p),
    length(beta), length(r_y), length(sd_x), length(sr2),
    length(ci_lower), length(ci_upper)
  ))
  pred_label <- function(i) {
    if (!is.null(predictor_names) && length(predictor_names) >= i)
      predictor_names[i] else as.character(i)
  }

  # Per-predictor dispatch
  for (i in seq_len(k_pred)) {
    pname <- pred_label(i)
    bi <- vget(b, i); sei <- vget(se, i); ti <- vget(t, i)
    pi_ <- vget(p, i); betai <- vget(beta, i)
    cili <- vget(ci_lower, i); ciui <- vget(ci_upper, i)
    ryi <- vget(r_y, i); sxi <- vget(sd_x, i); sri <- vget(sr2, i)

    if (has(bi) && has(sei)) {
      add(recalc_t_from_b_se(bi, sei, t = ti,
            b_digits = b_digits, se_digits = se_digits,
            t_digits = t_digits, rounding = rounding), pname)
    }
    if (has(ti) && has(df)) {
      add(recalc_p_from_t_df(ti, df, p = pi_, p_op = p_op,
            t_digits = t_digits, df_digits = 0,
            p_digits = p_digits, rounding = rounding), pname)
    }
    if (has(bi) && has(sei) && has(df)) {
      ci_arg <- if (has(cili) && has(ciui)) c(cili, ciui) else NULL
      add(recalc_ci_from_b_se(bi, sei, df, level = level, ci = ci_arg,
            b_digits = b_digits, se_digits = se_digits,
            df_digits = 0, ci_digits = ci_digits, rounding = rounding), pname)
    }
    if (has(bi) && has(sxi) && has(sd_y)) {
      add(recalc_beta_from_b(bi, sxi, sd_y, beta = betai,
            b_digits = b_digits, sd_x_digits = sd_x_digits,
            sd_y_digits = sd_y_digits, beta_digits = beta_digits, rounding = rounding), pname)
    }
    if (has(betai) && has(r2) && has(n) && has(k)) {
      add(recalc_t_bound_from_beta(betai, r2, n, k, t = ti,
            beta_digits = beta_digits, r2_digits = r2_digits,
            n_digits = n_digits, k_digits = k_digits,
            t_digits = t_digits, rounding = rounding), pname)
    }
    if (has(ti) && has(r2) && has(df)) {
      add(recalc_semipartial_r2_from_t(ti, r2, df, sr2 = sri,
            t_digits = t_digits, r2_digits = r2_digits,
            df_digits = 0, sr2_digits = sr2_digits, rounding = rounding), pname)
    }
  }

  # Model-level dispatch
  if (has(r2) && has(n) && has(k)) {
    add(recalc_f_from_r2(r2, n, k, f = f,
          r2_digits = r2_digits, n_digits = n_digits,
          k_digits = k_digits, f_digits = f_digits, rounding = rounding))
    add(recalc_adj_r2(r2, n, k, adj_r2 = adj_r2,
          r2_digits = r2_digits, n_digits = n_digits,
          k_digits = k_digits, adj_r2_digits = adj_r2_digits, rounding = rounding))
  }
  if (has(f) && has(n) && has(k)) {
    add(recalc_p_from_f(f, df1 = k, df2 = n - k - 1,
          p = p_model, p_op = p_model_op,
          f_digits = f_digits, df1_digits = k_digits,
          df2_digits = 0, p_digits = p_model_digits, rounding = rounding))
  } else if (has(r2) && has(n) && has(k)) {
    add(recalc_p_model_from_r2(r2, n, k, p_model = p_model,
          p_op = p_model_op,
          r2_digits = r2_digits, n_digits = n_digits,
          k_digits = k_digits, p_model_digits = p_model_digits, rounding = rounding))
  }
  if (!is.null(beta) && !is.null(r_y) && length(beta) == length(r_y)) {
    add(recalc_r2_from_betas_corrs(beta, r_y, r2 = r2,
          betas_digits = beta_digits, r_y_digits = r_y_digits,
          r2_digits = r2_digits, rounding = rounding))
  }
  if (!is.null(delta_r2_blocks)) {
    add(recalc_r2_from_blocks(delta_r2_blocks, r2 = r2,
          delta_r2_digits = delta_r2_blocks_digits,
          r2_digits = r2_digits, rounding = rounding))
  }

  # Cross-equation dispatch
  if (!is.null(beta) && length(beta) == 2 &&
      !is.null(r_y) && length(r_y) == 2) {
    add(recalc_r12_from_normal_eqs(beta[1], beta[2], r_y[1], r_y[2],
          beta1_digits = beta_digits, beta2_digits = beta_digits,
          r_y1_digits = r_y_digits, r_y2_digits = r_y_digits, rounding = rounding))
  }
  if (!is.null(r_y) && length(r_y) == 2 && has(r_12)) {
    add(recalc_betas_from_corrs_k2(r_y[1], r_y[2], r_12, betas = beta,
          r_y1_digits = r_y_digits, r_y2_digits = r_y_digits,
          r_12_digits = r_12_digits, betas_digits = beta_digits, rounding = rounding))
  }
  if (!is.null(r_y)) {
    add(recalc_r2_floor_from_correlations(r_y, r2 = r2,
          r_y_digits = r_y_digits, r2_digits = r2_digits, rounding = rounding))
  }
  if (!is.null(beta) && !is.null(r_y) && has(f) && has(n) && has(k) &&
      length(beta) == length(r_y)) {
    add(recalc_betar_from_f(beta, r_y, f, n, k,
          betas_digits = beta_digits, r_y_digits = r_y_digits,
          f_digits = f_digits, n_digits = n_digits, k_digits = k_digits, rounding = rounding))
  }

  # Misreporting diagnostics
  if (diagnostics) {
    # E1: beta vs b labeling
    if (!is.null(beta) && !is.null(r_y) && length(beta) == length(r_y) &&
        !is.null(sd_x) && length(sd_x) == length(beta) &&
        has(sd_y) && has(r2)) {
      add(diagnose_beta_label(
        coef_reported = beta, r_y = r_y,
        sd_x = sd_x, sd_y = sd_y, r2 = r2,
        coef_digits = beta_digits, r_y_digits = r_y_digits,
        sd_x_digits = sd_x_digits, sd_y_digits = sd_y_digits,
        r2_digits = r2_digits, rounding = rounding))
    }
    # E2: R^2 vs adjusted R^2 labeling
    if (has(r2) && !is.null(beta) && !is.null(r_y) &&
        length(beta) == length(r_y) && has(n) && has(k)) {
      add(diagnose_r2_label(
        r2_reported = r2, betas = beta, r_y = r_y, n = n, k = k,
        r2_digits = r2_digits, betas_digits = beta_digits,
        r_y_digits = r_y_digits, n_digits = n_digits, k_digits = k_digits, rounding = rounding))
    }
    # E3: one- vs two-tailed p, per predictor
    for (i in seq_len(k_pred)) {
      pname <- pred_label(i)
      ti <- vget(t, i); pi_ <- vget(p, i)
      if (has(ti) && has(df) && has(pi_)) {
        add(diagnose_p_tails(ti, df, p = pi_, p_op = p_op,
              t_digits = t_digits, df_digits = 0,
              p_digits = p_digits, rounding = rounding), pname)
      }
    }
  }

  dplyr::bind_rows(results)
}

#' Natural-language summary of an audit_regression() result
#'
#' Categorises rows into failed, passed, reconstruction-only, and diagnostic
#' groups; computes overall consistency verdict; and prints a human-readable
#' report listing each failure with reported value, recalculated interval, and
#' minimum interval gap. The overall conclusion accounts for diagnostic
#' resolutions: if a main check fails but exactly one labeling interpretation
#' (E1 / E2 / E3) is consistent with the report, the conclusion is
#' "inconsistent as reported, but resolved by [labeling]".
#'
#' @param audit Tibble returned by \code{audit_regression()}.
#' @param verbose Logical. If TRUE (default), \code{cat()}s the report.
#' @return Invisibly returns the report as a character scalar.
#' @examples
#' a <- audit_regression(
#'   n = 1923, k = 2, r2 = 0.41,
#'   beta = c(-0.39, 0.47), r_y = c(-0.53, 0.61),
#'   sd_y = 1.4, sd_x = c(1.2, 0.5),
#'   sd_y_digits = 1, sd_x_digits = 1
#' )
#' summarise_audit(a)
#' @export
summarise_audit <- function(audit, verbose = TRUE) {
  if (!"predictor" %in% names(audit)) audit$predictor <- "model"

  fmt_int <- function(lo, hi) {
    if (is.na(lo) || is.na(hi)) return("[NA, NA]")
    sprintf("[%s, %s]",
            formatC(lo, digits = 4, format = "g"),
            formatC(hi, digits = 4, format = "g"))
  }
  fmt_pred <- function(p) {
    if (is.na(p) || p == "model" || p == "") "" else sprintf(" [%s]", p)
  }
  fmt_num <- function(x) {
    if (is.na(x)) "NA" else formatC(x, digits = 4, format = "g")
  }
  interval_gap <- function(rep_lo, rep_hi, rec_lo, rec_hi) {
    max(0, rep_lo - rec_hi, rec_lo - rep_hi)
  }

  is_diag <- grepl("^E[123]", audit$check)
  diag <- audit[is_diag, , drop = FALSE]
  main <- audit[!is_diag, , drop = FALSE]

  has_cmp <- !is.na(main$consistent)
  passes <- main[has_cmp & main$consistent, , drop = FALSE]
  fails  <- main[has_cmp & !main$consistent, , drop = FALSE]
  recons <- main[!has_cmp, , drop = FALSE]

  # Group diagnostics by E# code and predictor
  diag_groups <- list()
  if (nrow(diag) > 0) {
    diag$prefix <- substr(diag$check, 1, 2)
    keys <- paste(diag$prefix, diag$predictor, sep = "::")
    for (key in unique(keys)) {
      g <- diag[keys == key, , drop = FALSE]
      cons <- g$consistent
      n_consistent <- sum(!is.na(cons) & cons)
      diag_groups[[length(diag_groups) + 1L]] <- list(
        prefix = g$prefix[1],
        predictor = g$predictor[1],
        resolved_by = if (n_consistent == 1L) g$check[!is.na(cons) & cons][1]
                      else NA_character_,
        all_inconsistent = all(!is.na(cons)) && !any(cons),
        group = g
      )
    }
  }
  any_resolved <- any(vapply(diag_groups,
    function(d) !is.na(d$resolved_by), logical(1)))

  L <- character()
  L <- c(L, "Regression report audit", strrep("=", 25), "")

  n_fail <- nrow(fails)
  n_cmp  <- n_fail + nrow(passes)

  if (n_fail == 0 && n_cmp > 0) {
    L <- c(L,
      "Conclusion: REPORTED QUANTITIES ARE MUTUALLY CONSISTENT within rounding.",
      sprintf("All %d direct check(s) pass; every recalculated interval overlaps the reported value's rounding interval.",
              n_cmp))
  } else if (n_fail == 0 && n_cmp == 0) {
    L <- c(L,
      "Conclusion: NO DIRECT CHECKS RUN.",
      "Insufficient inputs were supplied to compare any reported value against a recalculated quantity. See reconstructed quantities below.")
  } else if (any_resolved) {
    resolved_by <- vapply(diag_groups,
      function(d) if (!is.na(d$resolved_by)) sprintf("%s (%s)", d$prefix, d$resolved_by) else "",
      character(1))
    resolved_by <- resolved_by[resolved_by != ""]
    L <- c(L,
      "Conclusion: INCONSISTENT AS REPORTED, but resolved by labeling correction(s).",
      sprintf("%d of %d direct checks fail. The following labeling diagnostic(s) are uniquely consistent and would reconcile the report:",
              n_fail, n_cmp),
      paste0("  - ", resolved_by))
  } else {
    L <- c(L,
      "Conclusion: REPORTED QUANTITIES ARE NOT MUTUALLY CONSISTENT.",
      sprintf("%d of %d direct checks fail, and no standard label confusion (beta vs b, R^2 vs adjusted R^2, p tails) resolves the discrepancy.",
              n_fail, n_cmp))
  }
  L <- c(L, "")

  if (nrow(fails) > 0) {
    L <- c(L, "Failed checks:")
    for (i in seq_len(nrow(fails))) {
      r <- fails[i, ]
      gap <- interval_gap(r$reported_lower, r$reported_upper,
                          r$recalculated_lower, r$recalculated_upper)
      ref_line <- if (!is.na(r$reported)) {
        sprintf("    reported:    %s (rounding interval %s)",
                fmt_num(r$reported), fmt_int(r$reported_lower, r$reported_upper))
      } else {
        sprintf("    reference:   %s",
                fmt_int(r$reported_lower, r$reported_upper))
      }
      L <- c(L,
        sprintf("  %s%s", r$check, fmt_pred(r$predictor)),
        ref_line,
        sprintf("    recalculated:  %s",
                fmt_int(r$recalculated_lower, r$recalculated_upper)),
        sprintf("    interval gap: %s", fmt_num(gap)),
        ""
      )
    }
  }

  if (nrow(passes) > 0) {
    L <- c(L, "Passed checks:")
    for (i in seq_len(nrow(passes))) {
      r <- passes[i, ]
      rep_disp <- if (!is.na(r$reported)) fmt_num(r$reported)
                  else fmt_int(r$reported_lower, r$reported_upper)
      L <- c(L, sprintf(
        "  %s%s: reported %s, recalculated %s",
        r$check, fmt_pred(r$predictor),
        rep_disp,
        fmt_int(r$recalculated_lower, r$recalculated_upper)))
    }
    L <- c(L, "")
  }

  if (nrow(recons) > 0) {
    L <- c(L, "Reconstructed quantities (no reported value to compare):")
    for (i in seq_len(nrow(recons))) {
      r <- recons[i, ]
      L <- c(L, sprintf(
        "  %s%s: %s",
        r$check, fmt_pred(r$predictor),
        fmt_int(r$recalculated_lower, r$recalculated_upper)))
    }
    L <- c(L, "")
  }

  if (length(diag_groups) > 0) {
    L <- c(L, "Misreporting diagnostics:")
    for (d in diag_groups) {
      diag_name <- switch(d$prefix,
        E1 = "E1: standardized beta vs unstandardized b labeling",
        E2 = "E2: R^2 vs adjusted R^2 labeling",
        E3 = "E3: one- vs two-tailed p")
      L <- c(L, sprintf("  %s%s", diag_name, fmt_pred(d$predictor)))
      for (j in seq_len(nrow(d$group))) {
        r <- d$group[j, ]
        verdict <- if (is.na(r$consistent)) "n/a"
                   else if (r$consistent) "consistent"
                   else "not consistent"
        L <- c(L, sprintf(
          "    - %s: reported %s, recalculated %s -> %s",
          r$check,
          fmt_int(r$reported_lower, r$reported_upper),
          fmt_int(r$recalculated_lower, r$recalculated_upper),
          verdict))
      }
      verdict <- if (!is.na(d$resolved_by)) {
        sprintf("RESOLVES the discrepancy under '%s'.", d$resolved_by)
      } else if (d$all_inconsistent) {
        "does NOT explain the discrepancy (neither interpretation is consistent)."
      } else {
        "is uninformative (zero or multiple interpretations consistent)."
      }
      L <- c(L, sprintf("    Labeling verdict: %s", verdict), "")
    }
  }

  out <- paste(L, collapse = "\n")
  if (verbose) cat(out, "\n", sep = "")
  invisible(out)
}
