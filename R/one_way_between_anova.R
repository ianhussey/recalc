#' Multiverse of one-way between-subjects ANOVA from M/SD/N
#'
#' @description
#' Explore a multiverse of ANOVA results (F, p, partial eta^2, CIs) compatible
#' with reported group means, SDs, and Ns under different assumptions:
#' input rounding (±½ ULP), SE/SD confusion, CI methods, and output rounding.
#'
#' Implemented degrees of freedom:
#' \itemize{
#'   \item Input rounding: reported vs ± 0.5 ULP perturbation of means/SDs.
#'   \item SE/SD confusion: reported SDs vs interpreting SDs as SEs (SD = SE * sqrt(N)).
#'   \item CI method: none vs Wald-type CI on the correlation-equivalent scale.
#'   \item Output rounding: half-up / half-down / bankers / truncation
#'         (using \pkg{roundwork}).
#' }
#'
#' Not included by design:
#' \itemize{
#'   \item Alternative effect-size definitions (only partial eta^2).
#'   \item Unequal-variance ANOVA variants.
#'   \item Reconstruction paths other than computing F, p, and partial eta^2
#'         directly from M/SD/N.
#' }
#'
#' @param means Numeric vector of group means.
#' @param sds   Numeric vector of group SDs (or SEs in the SE/SD-confusion branch).
#' @param ns    Numeric vector of group sample sizes (integers).
#'
#' @param allow_input_rounding Logical vector of length 1 or 2. If TRUE is present,
#'   a branch using ± 0.5 ULP perturbations of means/SDs is created. If FALSE is
#'   present, a branch using the reported values as-is is created.
#'
#' @param se_sd_confusion Logical vector of length 1 or 2 indicating whether to
#'   include a branch that treats reported SDs as SEs (converted via SD = SE * sqrt(N)).
#'
#' @param ci_methods Character vector of CI methods for partial eta^2. Currently:
#'   \itemize{
#'     \item "none": no CIs (NA).
#'     \item "wald": Wald-type CI on the Fisher-z of the correlation equivalent.
#'   }
#'
#' @param output_rounding Character vector of rounding methods for outputs:
#'   \itemize{
#'     \item "half_up"
#'     \item "half_down"
#'     \item "bankers"
#'     \item "truncate"
#'   }
#'
#' @param alpha Significance level for CIs (default .05).
#'
#' @param digits_F   Number of digits for F when rounding outputs.
#' @param digits_p   Number of digits for p when rounding outputs.
#' @param digits_eta Number of digits for partial eta^2 when rounding outputs.
#' @param digits_ci  Number of digits for CI limits when rounding outputs.
#'
#' @return A list with two data frames:
#' \describe{
#'   \item{inputs}{Row per (scenario × group) with the "true" cell means, SDs, Ns,
#'                 and analytic assumptions (input rounding, SE/SD confusion).}
#'   \item{anova}{Row per (scenario × CI method × output-rounding) with F, df1, df2,
#'                p, partial eta^2, CIs (if any), and all analytic choices.}
#' }
#'
#' @importFrom stats pf qnorm
#' @importFrom roundwork round_up round_down round_bankers round_trunc
#' @export
multiverse_anova_oneway <- function(
    means,
    sds,
    ns,
    allow_input_rounding = c(FALSE, TRUE),
    se_sd_confusion      = c(FALSE, TRUE),
    ci_methods           = c("none", "wald"),
    output_rounding      = c("half_up", "half_down", "bankers", "truncate"),
    alpha                = 0.05,
    digits_F             = 2,
    digits_p             = 3,
    digits_eta           = 3,
    digits_ci            = 3
) {
  # ---- basic checks ----------------------------------------------------------
  if (length(means) != length(sds) || length(means) != length(ns)) {
    stop("means, sds, and ns must have the same length.")
  }
  if (any(ns <= 0)) {
    stop("All ns must be positive.")
  }
  if (length(means) < 2) {
    stop("At least two groups are required for one-way ANOVA.")
  }
  
  allow_input_rounding <- unique(allow_input_rounding)
  se_sd_confusion      <- unique(se_sd_confusion)
  ci_methods           <- unique(ci_methods)
  output_rounding      <- unique(output_rounding)
  
  K <- length(means)
  
  # ---- generate input scenarios: input rounding × SE/SD confusion -----------
  input_scenarios <- list()
  scenario_id     <- 0L
  
  for (ir in allow_input_rounding) {
    for (confuse in se_sd_confusion) {
      
      if (isTRUE(ir)) {
        # ± 0.5 ULP perturbations for all means and SDs
        mean_cands_list <- lapply(means, .candidate_values)
        sd_cands_list   <- lapply(sds,   .candidate_values)
        
        # Cartesian product
        means_grid <- expand.grid(mean_cands_list, KEEP.OUT.ATTRS = FALSE,
                                  stringsAsFactors = FALSE)
        colnames(means_grid) <- paste0("m", seq_len(K))
        
        sd_grid <- expand.grid(sd_cands_list, KEEP.OUT.ATTRS = FALSE,
                               stringsAsFactors = FALSE)
        colnames(sd_grid) <- paste0("sd", seq_len(K))
        
        grid <- cbind(means_grid, sd_grid)
        
      } else {
        # single scenario: reported means and SDs
        grid <- as.data.frame(
          c(
            as.list(setNames(as.list(means), paste0("m",  seq_len(K)))),
            as.list(setNames(as.list(sds),   paste0("sd", seq_len(K))))
          ),
          stringsAsFactors = FALSE
        )
      }
      
      n_scenarios <- nrow(grid)
      if (n_scenarios == 0) next
      
      # attach Ns
      for (j in seq_len(K)) {
        grid[[paste0("n", j)]] <- ns[j]
      }
      
      # SE/SD confusion
      if (isTRUE(confuse)) {
        for (j in seq_len(K)) {
          sd_col <- paste0("sd", j)
          n_col  <- paste0("n",  j)
          grid[[sd_col]] <- grid[[sd_col]] * sqrt(grid[[n_col]])
        }
      }
      
      # metadata
      grid$input_rounding  <- ir
      grid$se_sd_confusion <- confuse
      
      idx <- seq_len(n_scenarios)
      grid$scenario_id <- scenario_id + idx
      scenario_id      <- scenario_id + n_scenarios
      
      input_scenarios[[length(input_scenarios) + 1L]] <- grid
    }
  }
  
  if (!length(input_scenarios)) {
    stop("No input scenarios generated. Check arguments.")
  }
  
  inputs_df <- do.call(rbind, input_scenarios)
  
  # ---- long-format "inputs" data frame --------------------------------------
  inputs_long <- do.call(
    rbind,
    lapply(split(inputs_df, inputs_df$scenario_id), function(df) {
      sid   <- df$scenario_id[1]
      ir    <- df$input_rounding[1]
      conf  <- df$se_sd_confusion[1]
      row   <- df[1, , drop = FALSE]
      
      out_list <- vector("list", K)
      for (g in seq_len(K)) {
        out_list[[g]] <- data.frame(
          scenario_id      = sid,
          group            = g,
          mean             = row[[paste0("m",  g)]],
          sd               = row[[paste0("sd", g)]],
          n                = row[[paste0("n",  g)]],
          input_rounding   = ir,
          se_sd_confusion  = conf,
          stringsAsFactors = FALSE
        )
      }
      do.call(rbind, out_list)
    })
  )
  
  # ---- ANOVA calculations per scenario × CI method × rounding ---------------
  anova_rows <- list()
  
  for (i in seq_len(nrow(inputs_df))) {
    scen <- inputs_df[i, ]
    
    m_vec <- as.numeric(scen[ , paste0("m",  seq_len(K))])
    sd_vec <- as.numeric(scen[ , paste0("sd", seq_len(K))])
    n_vec  <- as.numeric(scen[ , paste0("n",  seq_len(K))])
    
    N_total    <- sum(n_vec)
    grand_mean <- sum(n_vec * m_vec) / N_total
    
    SS_effect <- sum(n_vec * (m_vec - grand_mean)^2)
    SS_error  <- sum((n_vec - 1) * (sd_vec^2))
    
    df1 <- K - 1L
    df2 <- N_total - K
    
    if (SS_error <= 0 || df2 <= 0) {
      F_val   <- NA_real_
      p_val   <- NA_real_
      eta_p2  <- NA_real_
    } else {
      MS_effect <- SS_effect / df1
      MS_error  <- SS_error / df2
      
      F_val  <- MS_effect / MS_error
      p_val  <- 1 - stats::pf(F_val, df1, df2)
      eta_p2 <- SS_effect / (SS_effect + SS_error)
    }
    
    for (ci_method in ci_methods) {
      
      if (is.na(eta_p2) || ci_method == "none") {
        ci_low  <- NA_real_
        ci_high <- NA_real_
      } else if (ci_method == "wald") {
        ci      <- .ci_eta_wald(eta_p2, df1, df2, alpha)
        ci_low  <- ci[1]
        ci_high <- ci[2]
      } else {
        stop("Unknown ci_method: ", ci_method)
      }
      
      for (orule in output_rounding) {
        
        F_out     <- .round_output_rw(F_val,    digits_F,   orule)
        p_out     <- .round_output_rw(p_val,    digits_p,   orule)
        eta_out   <- .round_output_rw(eta_p2,   digits_eta, orule)
        ci_low_o  <- .round_output_rw(ci_low,   digits_ci,  orule)
        ci_high_o <- .round_output_rw(ci_high,  digits_ci,  orule)
        
        anova_rows[[length(anova_rows) + 1L]] <- data.frame(
          scenario_id      = scen$scenario_id,
          input_rounding   = scen$input_rounding,
          se_sd_confusion  = scen$se_sd_confusion,
          ci_method        = ci_method,
          output_rounding  = orule,
          df1              = df1,
          df2              = df2,
          F                = F_out,
          p                = p_out,
          eta2_partial     = eta_out,
          eta2_ci_low      = ci_low_o,
          eta2_ci_high     = ci_high_o,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  anova_df <- do.call(rbind, anova_rows)
  
  # ---- return ---------------------------------------------------------------
  list(
    inputs = inputs_long,
    anova  = anova_df
  )
}

#' Internal: count decimal places in a numeric vector
#'
#' @keywords internal
.decimal_places <- function(x) {
  sapply(x, function(xx) {
    if (is.na(xx)) return(NA_integer_)
    ch <- format(xx, scientific = FALSE, trim = TRUE)
    if (!grepl("\\.", ch, fixed = TRUE)) {
      0L
    } else {
      nchar(strsplit(ch, ".", fixed = TRUE)[[1]][2])
    }
  })
}

#' Internal: ± 0.5 ULP candidates for one reported value
#'
#' @keywords internal
.candidate_values <- function(x) {
  d   <- .decimal_places(x)
  ulp <- 0.5 * 10^(-d)
  c(x - ulp, x, x + ulp)
}

#' Internal: rounding wrapper using {roundwork}
#'
#' @param x numeric
#' @param digits integer
#' @param method one of "half_up", "half_down", "bankers", "truncate"
#'
#' @keywords internal
.round_output_rw <- function(x, digits, method = c("half_up", "half_down",
                                                   "bankers", "truncate")) {
  method <- match.arg(method)
  
  if (all(is.na(x)))    return(x)
  
  # assume these functions exist in {roundwork}; adjust names if needed
  fun <- switch(
    method,
    half_up   = roundwork::round_up,
    half_down = roundwork::round_down,
    bankers   = round,
    truncate  = roundwork::round_trunc,  
    stop("Unknown rounding method: ", method)
  )
  
  fun(x, digits = digits)
}

#' Internal: Wald-type CI for partial eta^2 via correlation equivalence
#'
#' @keywords internal
.ci_eta_wald <- function(eta2, df1, df2, alpha = 0.05) {
  if (is.na(eta2) || eta2 <= 0 || eta2 >= 1) {
    return(c(NA_real_, NA_real_))
  }
  
  r <- sqrt(eta2)
  # Effective N for the correlation equivalent
  N_eff <- df1 + df2 + 1L
  if (N_eff <= 3) {
    return(c(NA_real_, NA_real_))
  }
  
  z      <- 0.5 * log((1 + r) / (1 - r))
  se_z   <- 1 / sqrt(N_eff - 3)
  z_crit <- stats::qnorm(1 - alpha / 2)
  
  z_low  <- z - z_crit * se_z
  z_high <- z + z_crit * se_z
  
  r_low  <- (exp(2 * z_low) - 1) / (exp(2 * z_low) + 1)
  r_high <- (exp(2 * z_high) - 1) / (exp(2 * z_high) + 1)
  
  c(r_low^2, r_high^2)
}