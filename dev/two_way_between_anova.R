#' Multiverse of two-way between-subjects ANOVA from M/SD/N
#'
#' @param means Matrix a × b of cell means.
#' @param sds   Matrix a × b of cell SDs (or SEs under confusion).
#' @param ns    Matrix a × b of cell sample sizes.
#'
#' @details
#' Implements fully crossed multiverse:
#'   input rounding × SE/SD confusion × CI method × output rounding.
#'
#' Effects returned separately:
#'   * main effect A
#'   * main effect B
#'   * interaction A×B
#'
#' Only uses:  
#'   - partial η²  
#'   - Wald CIs (safe from summary stats)  
#'
#' @return List with: inputs, anova_A, anova_B, anova_AB.
#'
#' @export
two_way_between_anova <- function(
    means,
    sds,
    ns,
    allow_input_rounding = c(FALSE, TRUE),
    se_sd_confusion      = c(FALSE, TRUE),
    ci_methods           = c("none", "wald"),
    output_rounding      = c("none", "half_up", "half_down", "bankers", "truncate"),
    alpha                = 0.05,
    digits_F             = 2,
    digits_p             = 3,
    digits_eta           = 3,
    digits_ci            = 3
) {
  
  # ---------------------------------------------------------------------------
  # checks
  # ---------------------------------------------------------------------------
  if (!is.matrix(means) || !is.matrix(sds) || !is.matrix(ns))
    stop("means, sds, ns must all be matrices with identical dimensions.")
  
  if (!all(dim(means) == dim(sds)) || !all(dim(means) == dim(ns)))
    stop("means, sds, ns must have identical dimensions (a × b).")
  
  a <- nrow(means)
  b <- ncol(means)
  
  # flatten cell indices for easier manipulation
  idx <- expand.grid(A = seq_len(a), B = seq_len(b))
  
  allow_input_rounding <- unique(allow_input_rounding)
  se_sd_confusion      <- unique(se_sd_confusion)
  ci_methods           <- unique(ci_methods)
  output_rounding      <- unique(output_rounding)
  
  # ---------------------------------------------------------------------------
  # generate input scenarios
  # ---------------------------------------------------------------------------
  input_scenarios <- list()
  scenario_id     <- 0L
  
  for (ir in allow_input_rounding) {
    for (confuse in se_sd_confusion) {
      
      if (isTRUE(ir)) {
        mean_cand_list <- lapply(as.vector(means), .candidate_values)
        sd_cand_list   <- lapply(as.vector(sds),   .candidate_values)
        
        mean_grid <- expand.grid(mean_cand_list, KEEP.OUT.ATTRS = FALSE)
        sd_grid   <- expand.grid(sd_cand_list,   KEEP.OUT.ATTRS = FALSE)
        
        colnames(mean_grid) <- paste0("m", seq_len(nrow(idx)))
        colnames(sd_grid)   <- paste0("sd", seq_len(nrow(idx)))
        
        grid <- cbind(mean_grid, sd_grid)
      } else {
        mvec <- as.vector(means)
        sdv  <- as.vector(sds)
        grid <- as.data.frame(
          c(
            setNames(as.list(mvec), paste0("m", seq_len(nrow(idx)))),
            setNames(as.list(sdv),  paste0("sd", seq_len(nrow(idx))))
          )
        )
      }
      
      # attach ns
      for (j in seq_len(nrow(idx))) {
        grid[[paste0("n", j)]] <- as.vector(ns)[j]
      }
      
      # SE/SD confusion
      if (isTRUE(confuse)) {
        for (j in seq_len(nrow(idx))) {
          sdcol <- paste0("sd", j)
          ncol  <- paste0("n",  j)
          grid[[sdcol]] <- grid[[sdcol]] * sqrt(grid[[ncol]])
        }
      }
      
      n_scen <- nrow(grid)
      grid$scenario_id     <- seq(from = scenario_id + 1L, length.out = n_scen)
      grid$input_rounding  <- ir
      grid$se_sd_confusion <- confuse
      scenario_id          <- scenario_id + n_scen
      
      input_scenarios[[length(input_scenarios) + 1L]] <- grid
    }
  }
  
  inputs_df <- do.call(rbind, input_scenarios)
  
  # long-format inputs
  inputs_long <- do.call(
    rbind,
    lapply(split(inputs_df, inputs_df$scenario_id), function(df) {
      sid  <- df$scenario_id[1]
      ir   <- df$input_rounding[1]
      conf <- df$se_sd_confusion[1]
      row  <- df[1, , drop = FALSE]
      
      out <- vector("list", nrow(idx))
      for (j in seq_len(nrow(idx))) {
        out[[j]] <- data.frame(
          scenario_id     = sid,
          A               = idx$A[j],
          B               = idx$B[j],
          mean            = row[[paste0("m", j)]],
          sd              = row[[paste0("sd", j)]],
          n               = row[[paste0("n", j)]],
          input_rounding  = ir,
          se_sd_confusion = conf
        )
      }
      do.call(rbind, out)
    })
  )
  
  # ---------------------------------------------------------------------------
  # Function to compute SS for A, B, AB for one scenario
  # ---------------------------------------------------------------------------
  compute_two_way_SS <- function(mvec, sdvec, nvec) {
    # reshape into a×b
    M  <- matrix(mvec, a, b, byrow = FALSE)
    SD <- matrix(sdvec, a, b, byrow = FALSE)
    N  <- matrix(nvec,  a, b, byrow = FALSE)
    
    N_total <- sum(N)
    grand_mean <- sum(N * M) / N_total
    
    # marginal means
    M_A <- rowSums(N * M) / rowSums(N)
    M_B <- colSums(N * M) / colSums(N)
    
    # SS_effect_A
    SS_A <- sum(rowSums(N) * (M_A - grand_mean)^2)
    
    # SS_effect_B
    SS_B <- sum(colSums(N) * (M_B - grand_mean)^2)
    
    # SS_effect_AB
    SS_cells <- sum(N * (M - grand_mean)^2)
    SS_AB <- SS_cells - SS_A - SS_B
    
    # SS_error
    SS_error <- sum((N - 1) * (SD^2))
    
    list(
      SS_A = SS_A,
      SS_B = SS_B,
      SS_AB = SS_AB,
      SS_error = SS_error,
      N_total = N_total
    )
  }
  
  # ---------------------------------------------------------------------------
  # ANALYSE ALL SCENARIOS
  # ---------------------------------------------------------------------------
  out_A  <- list()
  out_B  <- list()
  out_AB <- list()
  
  for (i in seq_len(nrow(inputs_df))) {
    scen <- inputs_df[i, ]
    
    mvec <- as.numeric(scen[paste0("m", seq_len(nrow(idx)))])
    sdv  <- as.numeric(scen[paste0("sd", seq_len(nrow(idx)))])
    nv   <- as.numeric(scen[paste0("n", seq_len(nrow(idx)))])
    
    SS <- compute_two_way_SS(mvec, sdv, nv)
    
    dfA  <- a - 1L
    dfB  <- b - 1L
    dfAB <- (a - 1L) * (b - 1L)
    dfE  <- SS$N_total - a*b
    
    effects <- list(
      A  = list(SS = SS$SS_A,  df = dfA),
      B  = list(SS = SS$SS_B,  df = dfB),
      AB = list(SS = SS$SS_AB, df = dfAB)
    )
    
    SS_error <- SS$SS_error
    
    for (ci_method in ci_methods) {
      for (orule in output_rounding) {
        
        for (effect_name in c("A","B","AB")) {
          
          SSx <- effects[[effect_name]]$SS
          dfx <- effects[[effect_name]]$df
          
          if (SS_error <= 0 || dfE <= 0) {
            Fval  <- NA_real_
            pval  <- NA_real_
            etax  <- NA_real_
            ci_l  <- NA_real_
            ci_h  <- NA_real_
          } else {
            MSx <- SSx / dfx
            MSe <- SS_error / dfE
            Fval <- MSx / MSe
            pval <- 1 - stats::pf(Fval, dfx, dfE)
            etax <- SSx / (SSx + SS_error)
            
            if (ci_method == "wald") {
              ci <- .ci_eta_wald(etax, dfx, dfE, alpha)
              ci_l <- ci[1]
              ci_h <- ci[2]
            } else {
              ci_l <- ci_h <- NA_real_
            }
          }
          
          df_out <- data.frame(
            scenario_id     = scen$scenario_id,
            effect          = effect_name,
            input_rounding  = scen$input_rounding,
            se_sd_confusion = scen$se_sd_confusion,
            ci_method       = ci_method,
            output_rounding = orule,
            df1             = dfx,
            df2             = dfE,
            F               = .round_output_rw(Fval,  digits_F, orule),
            p               = .round_output_rw(pval,  digits_p, orule),
            eta2_partial    = .round_output_rw(etax, digits_eta, orule),
            eta2_ci_low     = .round_output_rw(ci_l, digits_ci, orule),
            eta2_ci_high    = .round_output_rw(ci_h, digits_ci, orule)
          )
          
          if (effect_name == "A")  out_A[[length(out_A)+1]]   <- df_out
          if (effect_name == "B")  out_B[[length(out_B)+1]]   <- df_out
          if (effect_name == "AB") out_AB[[length(out_AB)+1]] <- df_out
        }
      }
    }
  }
  
  list(
    inputs  = inputs_long,
    anova_A = dplyr::bind_rows(out_A),
    anova_B = dplyr::bind_rows(out_B),
    anova_AB = dplyr::bind_rows(out_AB)
  )
}


#' @keywords internal
.decimal_places <- function(x) {
  sapply(x, function(xx) {
    if (is.na(xx)) return(NA_integer_)
    ch <- format(xx, scientific = FALSE, trim = TRUE)
    if (!grepl("\\.", ch, fixed = TRUE)) 0L else {
      nchar(strsplit(ch, ".", fixed = TRUE)[[1]][2])
    }
  })
}

#' @keywords internal
.candidate_values <- function(x) {
  d   <- .decimal_places(x)
  ulp <- 0.5 * 10^(-d)
  c(x - ulp, x, x + ulp)
}

#' @keywords internal
.round_output_rw <- function(x, digits, method) {
  method <- match.arg(method, c("none","half_up","half_down","bankers","truncate"))
  if (method == "none") return(x)
  if (all(is.na(x)))    return(x)
  
  fun <- switch(
    method,
    half_up   = roundwork::round_up,
    half_down = roundwork::round_down,
    bankers   = round,
    truncate  = roundwork::round_trunc
  )
  fun(x, digits = digits)
}

#' @keywords internal
.ci_eta_wald <- function(eta2, df1, df2, alpha = 0.05) {
  if (is.na(eta2) || eta2 <= 0 || eta2 >= 1) return(c(NA_real_, NA_real_))
  r <- sqrt(eta2)
  N_eff <- df1 + df2 + 1L
  if (N_eff <= 3) return(c(NA_real_, NA_real_))
  
  z      <- 0.5 * log((1 + r)/(1 - r))
  zcrit  <- stats::qnorm(1 - alpha/2)
  se_z   <- 1/sqrt(N_eff - 3)
  
  z_low  <- z - zcrit * se_z
  z_high <- z + zcrit * se_z
  
  r_low  <- (exp(2*z_low) - 1)/(exp(2*z_low) + 1)
  r_high <- (exp(2*z_high) - 1)/(exp(2*z_high) + 1)
  
  c(r_low^2, r_high^2)
}

