
# todo
# next add t_test_from_summary_naudet.r or t-test.R to also calculate t test results and p values from summary stats, with a range of p values for rounding, welch/student, etc.

identify_d_method <- function(m1, m2, sd1, sd2, n1, n2,
                              t = NULL, df = NULL,
                              d_est, d_ci_lower, d_ci_upper,
                              d_digits = 2,
                              alpha = 0.05) {
  require(roundwork)
  
  # Coerce reported values to numeric if they came in as characters
  d_est      <- as.numeric(d_est)
  d_ci_lower <- as.numeric(d_ci_lower)
  d_ci_upper <- as.numeric(d_ci_upper)
  
  # Helper: infer how many decimal places the user supplied
  get_digits <- function(x) {
    sx <- sub("^-", "", as.character(x))
    if (!grepl("\\.", sx)) return(0L)
    nchar(sub("^[^.]*\\.", "", sx))
  }
  
  # Helper: apply a +/- half-ulp adjustment to a value
  adjust_value <- function(x, step, code) {
    if (step == 0) return(x)
    if (code == "minus") return(x - step)
    if (code == "plus")  return(x + step)
    x  # "reported"
  }
  
  # Noncentral t CI for the noncentrality parameter delta
  nct_ci <- function(t_obs, df, alpha = 0.05, max_ncp = 1000) {
    # Returns c(delta_L, delta_U); NA if root finding fails
    fL <- function(delta) stats::pt(t_obs, df = df, ncp = delta) - alpha / 2
    fU <- function(delta) stats::pt(t_obs, df = df, ncp = delta) - (1 - alpha / 2)
    
    lower <- -max_ncp
    upper <-  max_ncp
    
    # lower limit
    fL_low  <- fL(lower)
    fL_high <- fL(upper)
    if (is.na(fL_low) || is.na(fL_high) || fL_low * fL_high > 0) {
      return(c(NA_real_, NA_real_))
    }
    delta_L <- uniroot(fL, lower = lower, upper = upper)$root
    
    # upper limit
    fU_low  <- fU(lower)
    fU_high <- fU(upper)
    if (is.na(fU_low) || is.na(fU_high) || fU_low * fU_high > 0) {
      return(c(NA_real_, NA_real_))
    }
    delta_U <- uniroot(fU, lower = lower, upper = upper)$root
    
    c(delta_L, delta_U)
  }
  
  # What information do we have?
  have_summary <- all(!is.null(c(m1, m2, sd1, sd2, n1, n2))) &&
    all(!is.na   (c(m1, m2, sd1, sd2, n1, n2)))
  
  have_t <- !is.null(t) && !is.null(df) &&
    !is.na(t)   && !is.na(df)   &&
    all(!is.null(c(n1, n2))) &&
    all(!is.na   (c(n1, n2)))
  
  if (!have_summary && !have_t) {
    stop("Provide either M/SD/N (m1, m2, sd1, sd2, n1, n2) and/or t and df.")
  }
  
  results <- list()
  idx <- 1
  
  ## 0) Reported values row ----
  results[[idx]] <- data.frame(
    source              = "reported",
    direction           = NA_character_,
    es_type             = NA_character_,
    ci_method           = NA_character_,
    rounding            = NA_character_,
    input_adj_stats     = NA_character_,
    input_adj_tdf       = NA_character_,
    m1_used             = NA_real_,
    m2_used             = NA_real_,
    sd1_used            = NA_real_,
    sd2_used            = NA_real_,
    es_unrounded        = NA_real_,
    ci_lower_unrounded  = NA_real_,
    ci_upper_unrounded  = NA_real_,
    est_rounded         = d_est,
    ci_lower_rounded    = d_ci_lower,
    ci_upper_rounded    = d_ci_upper,
    match_est           = NA,
    match_ci_lower      = NA,
    match_ci_upper      = NA,
    match_all           = NA,
    stringsAsFactors    = FALSE
  )
  idx <- idx + 1
  
  ci_methods <- c("wald_z", "wald_t", "welch_t", "welch_z", "nct")
  
  #### 1) From summary statistics (means, SDs, Ns) ####
  if (have_summary) {
    
    # Step sizes = half of one unit in last reported decimal place
    dig_m1  <- get_digits(m1)
    dig_m2  <- get_digits(m2)
    dig_sd1 <- get_digits(sd1)
    dig_sd2 <- get_digits(sd2)
    
    step_m1  <- 0.5 * 10^(-dig_m1)
    step_m2  <- 0.5 * 10^(-dig_m2)
    step_sd1 <- 0.5 * 10^(-dig_sd1)
    step_sd2 <- 0.5 * 10^(-dig_sd2)
    
    # Codes for input rounding assumptions applied uniformly to all stats
    adj_codes <- c("reported", "minus", "plus")
    
    for (adj_stats in adj_codes) {
      
      m1_star  <- adjust_value(m1,  step_m1,  adj_stats)
      m2_star  <- adjust_value(m2,  step_m2,  adj_stats)
      sd1_star <- adjust_value(sd1, step_sd1, adj_stats)
      sd2_star <- adjust_value(sd2, step_sd2, adj_stats)
      
      # Skip impossible SDs
      if (sd1_star <= 0 || sd2_star <= 0) next
      
      for (direction in c("m1_minus_m2", "m2_minus_m1")) {
        
        diff_mean <- if (direction == "m1_minus_m2") {
          m1_star - m2_star
        } else {
          m2_star - m1_star
        }
        
        # Pooled SD and basic quantities
        df_s  <- n1 + n2 - 2
        N     <- n1 + n2
        sp    <- sqrt(((n1 - 1) * sd1_star^2 + (n2 - 1) * sd2_star^2) / df_s)
        if (!is.finite(sp) || sp <= 0) next
        
        d_raw <- diff_mean / sp
        
        # Welch df (for unequal variances)
        var1 <- sd1_star^2
        var2 <- sd2_star^2
        num_w <- (var1 / n1 + var2 / n2)^2
        den_w <- (var1^2 / (n1^2 * (n1 - 1))) +
          (var2^2 / (n2^2 * (n2 - 1)))
        df_w <- num_w / den_w
        if (!is.finite(df_w) || df_w <= 0) df_w <- NA_real_
        
        # t statistic using pooled SD (needed for nct)
        t_pooled <- diff_mean / (sp * sqrt(1 / n1 + 1 / n2))
        
        # Hedges' J small-sample correction
        J_s   <- 1 - 3 / (4 * df_s - 1)
        g_raw <- J_s * d_raw
        
        # Candidate effect sizes: Cohen's d vs Hedges' g
        for (es_type in c("d", "g")) {
          
          es <- if (es_type == "d") d_raw else g_raw
          
          # SE under pooled df and Welch df
          se_pooled <- sqrt(N / (n1 * n2) + es^2 / (2 * df_s))
          se_welch  <- if (!is.na(df_w)) {
            sqrt(N / (n1 * n2) + es^2 / (2 * df_w))
          } else {
            NA_real_
          }
          
          for (ci_method in ci_methods) {
            
            if (ci_method == "nct") {
              # Noncentral t CI based on pooled t and df_s
              delta_ci <- nct_ci(t_pooled, df_s, alpha = alpha)
              if (any(is.na(delta_ci))) next
              
              # Transform delta to d; then to g if needed
              fac_d <- sqrt(1 / n1 + 1 / n2)
              dL_raw <- delta_ci[1] * fac_d
              dU_raw <- delta_ci[2] * fac_d
              
              if (es_type == "d") {
                ci_lower <- dL_raw
                ci_upper <- dU_raw
              } else {
                ci_lower <- J_s * dL_raw
                ci_upper <- J_s * dU_raw
              }
              
            } else {
              # Wald-type CIs
              if (ci_method %in% c("wald_z", "wald_t")) {
                se_use   <- se_pooled
                df_for_t <- df_s
              } else { # welch_t, welch_z
                se_use   <- if (!is.na(se_welch)) se_welch else se_pooled
                df_for_t <- if (!is.na(df_w)) df_w else df_s
              }
              
              if (ci_method %in% c("wald_z", "welch_z")) {
                crit <- stats::qnorm(1 - alpha / 2)
              } else {
                crit <- stats::qt(1 - alpha / 2, df = df_for_t)
              }
              
              ci_lower <- es - crit * se_use
              ci_upper <- es + crit * se_use
            }
            
            for (rounding in c("half_up", "half_down")) {
              
              round_fun <- if (rounding == "half_up") {
                function(x) round_up(x, d_digits)
              } else {
                function(x) round_down(x, d_digits)
              }
              
              est_r   <- round_fun(es)
              lower_r <- round_fun(ci_lower)
              upper_r <- round_fun(ci_upper)
              
              results[[idx]] <- data.frame(
                source              = "summary",
                direction           = direction,
                es_type             = es_type,
                ci_method           = ci_method,
                rounding            = rounding,
                input_adj_stats     = adj_stats,
                input_adj_tdf       = NA_character_,
                m1_used             = m1_star,
                m2_used             = m2_star,
                sd1_used            = sd1_star,
                sd2_used            = sd2_star,
                es_unrounded        = es,
                ci_lower_unrounded  = ci_lower,
                ci_upper_unrounded  = ci_upper,
                est_rounded         = est_r,
                ci_lower_rounded    = lower_r,
                ci_upper_rounded    = upper_r,
                match_est           = isTRUE(all.equal(est_r, d_est)),
                match_ci_lower      = isTRUE(all.equal(lower_r, d_ci_lower)),
                match_ci_upper      = isTRUE(all.equal(upper_r, d_ci_upper)),
                match_all           = (est_r == d_est &&
                                         lower_r == d_ci_lower &&
                                         upper_r == d_ci_upper),
                stringsAsFactors    = FALSE
              )
              
              idx <- idx + 1
            }
          }
        }
      }
    }
  }
  
  #### 2) From t statistic and df ####
  if (have_t) {
    t_val <- as.numeric(t)
    df_t  <- as.numeric(df)
    N     <- n1 + n2
    df_s  <- n1 + n2 - 2  # pooled df
    
    # Step sizes for t and df (half-ulp of reported precision)
    dig_t  <- get_digits(t_val)
    dig_df <- get_digits(df_t)
    
    step_t  <- 0.5 * 10^(-dig_t)
    step_df <- 0.5 * 10^(-dig_df)
    
    adj_codes <- c("reported", "minus", "plus")
    
    for (adj_tdf in adj_codes) {
      
      t_star  <- adjust_value(t_val,  step_t,  adj_tdf)
      df_star <- adjust_value(df_t,  step_df, adj_tdf)
      if (!is.finite(df_star) || df_star <= 0) next
      
      # Classical relationship for independent groups with pooled SD:
      # d = t * sqrt(1/n1 + 1/n2)
      fac_d   <- sqrt(1 / n1 + 1 / n2)
      d_from_t <- t_star * fac_d
      
      # Allow possible sign flip (authors may have reversed group labels)
      for (sign_flip in c(1, -1)) {
        
        es_base  <- d_from_t * sign_flip
        direction <- if (sign_flip == 1) "t_sign" else "neg_t_sign"
        
        # Hedges' J for this df
        J_t   <- 1 - 3 / (4 * df_star - 1)
        g_base <- J_t * es_base
        
        for (es_type in c("d", "g")) {
          
          es <- if (es_type == "d") es_base else g_base
          
          # SE under pooled df_s and reported df_star
          se_pooled <- sqrt(N / (n1 * n2) + es^2 / (2 * df_s))
          se_welch  <- sqrt(N / (n1 * n2) + es^2 / (2 * df_star))
          
          for (ci_method in ci_methods) {
            
            if (ci_method == "nct") {
              # Noncentral t CI based on t_star and df_star
              delta_ci <- nct_ci(t_star, df_star, alpha = alpha)
              if (any(is.na(delta_ci))) next
              
              dL_raw <- delta_ci[1] * fac_d
              dU_raw <- delta_ci[2] * fac_d
              
              if (es_type == "d") {
                ci_lower <- dL_raw
                ci_upper <- dU_raw
              } else {
                ci_lower <- J_t * dL_raw
                ci_upper <- J_t * dU_raw
              }
              
            } else {
              if (ci_method %in% c("wald_z", "wald_t")) {
                se_use   <- se_pooled
                df_for_t <- df_s
              } else { # welch_t, welch_z
                se_use   <- se_welch
                df_for_t <- df_star
              }
              
              if (ci_method %in% c("wald_z", "welch_z")) {
                crit <- stats::qnorm(1 - alpha / 2)
              } else {
                crit <- stats::qt(1 - alpha / 2, df = df_for_t)
              }
              
              ci_lower <- es - crit * se_use
              ci_upper <- es + crit * se_use
            }
            
            for (rounding in c("half_up", "half_down")) {
              
              round_fun <- if (rounding == "half_up") {
                function(x) round_up(x, d_digits)
              } else {
                function(x) round_down(x, d_digits)
              }
              
              est_r   <- round_fun(es)
              lower_r <- round_fun(ci_lower)
              upper_r <- round_fun(ci_upper)
              
              results[[idx]] <- data.frame(
                source              = "t_df",
                direction           = direction,
                es_type             = es_type,
                ci_method           = ci_method,
                rounding            = rounding,
                input_adj_stats     = NA_character_,
                input_adj_tdf       = adj_tdf,
                m1_used             = NA_real_,
                m2_used             = NA_real_,
                sd1_used            = NA_real_,
                sd2_used            = NA_real_,
                es_unrounded        = es,
                ci_lower_unrounded  = ci_lower,
                ci_upper_unrounded  = ci_upper,
                est_rounded         = est_r,
                ci_lower_rounded    = lower_r,
                ci_upper_rounded    = upper_r,
                match_est           = isTRUE(all.equal(est_r, d_est)),
                match_ci_lower      = isTRUE(all.equal(lower_r, d_ci_lower)),
                match_ci_upper      = isTRUE(all.equal(upper_r, d_ci_upper)),
                match_all           = (est_r == d_est &&
                                         lower_r == d_ci_lower &&
                                         upper_r == d_ci_upper),
                stringsAsFactors    = FALSE
              )
              
              idx <- idx + 1
            }
          }
        }
      }
    }
  }
  
  out <- do.call(rbind, results)
  
  # Put exact matches (if any) at the top (reported row has match_all = NA, so is unaffected)
  out <- out[order(!out$match_all,
                   !out$match_est,
                   !out$match_ci_lower,
                   !out$match_ci_upper), ]
  
  list(
    reproduced = any(out$match_all, na.rm = TRUE),
    results    = out
  )
}


# example
res <- identify_d_method(
  m1 = 10.3, 
  m2 = 8.7,
  sd1 = 3.1, 
  sd2 = 2.8,
  n1 = 50, 
  n2 = 48,
  t = 1.34, 
  df = 97,
  d_est = 0.55,
  d_ci_lower = 0.20,
  d_ci_upper = 0.90,
  d_digits = 2
)$results |>
  as_tibble() |>
  mutate(original = source=="reported") |>
  arrange(est_rounded, original) |>
  rownames_to_column(var = "rowname") |>
  mutate(rowname = as.numeric(as.character(rowname)))
  
res |>
  filter(source != "reported" & direction %in% c("m1_minus_m2", "t_sign")) |>
  ggplot(aes(y = rowname, 
             x = est_rounded, xmin = ci_lower_rounded, xmax = ci_upper_rounded,
             color = original)) +
  ggstance::geom_linerangeh() +
  geom_point()


res <- identify_d_method(
  m1 = 10.3, 
  m2 = 8.7,
  sd1 = 3.1, 
  sd2 = 2.8,
  n1 = 50, 
  n2 = 48,
  #t = 1.34, 
  #df = 97,
  d_est = 0.55,
  d_ci_lower = 0.20,
  d_ci_upper = 0.90,
  d_digits = 2
)$results |>
  as_tibble() |>
  filter(source != "reported" & direction %in% c("m1_minus_m2", "t_sign")) 

res |>
  summarize(min_est_rounded = min(est_rounded),
            max_est_rounded = max(est_rounded))

