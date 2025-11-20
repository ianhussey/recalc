
# choices for p values and Cohen's d are independent, given people can't be assumed to be consistent
# add summary df with min/max for cohens d and its CIs and p values. how to handle cis with d ties? eg widest or narrowest?

# requires you to specify Ns and M+SD and/or t/df. Could write it to calculate d+p assuming even split, but this is more a meta-analysis method of making a reasonable estimation of effect size with available information rather than recalculate p/d.
# a seperate function could ask the question what allocation ratio allows you to recalculate p/d given rounding etc. to assess whether that allocation ratio is plausible. 

# maybe t/df and M/SD approaches should be separated into different functions. potential lack of agreement between them is hidden here. 

identify_d_method <- function(n1 = NULL, n2 = NULL,
                              m1 = NULL, m2 = NULL, sd1 = NULL, sd2 = NULL, 
                              t = NULL, df = NULL,
                              d_est = NULL, d_ci_lower = NULL, d_ci_upper = NULL,
                              d_digits = 2,
                              p_est = NULL, p_digits = 3,
                              alpha = 0.05) {
  require(roundwork)
  
  ## ------------------------------------------------------------------
  ## 0) Coerce reported values and set up helpers
  ## ------------------------------------------------------------------
  d_est_num      <- if (!is.null(d_est))      as.numeric(d_est)      else NA_real_
  d_ci_lower_num <- if (!is.null(d_ci_lower)) as.numeric(d_ci_lower) else NA_real_
  d_ci_upper_num <- if (!is.null(d_ci_upper)) as.numeric(d_ci_upper) else NA_real_
  p_est_num      <- if (!is.null(p_est))      as.numeric(p_est)      else NA_real_
  
  # Helper: infer how many decimal places the user supplied
  get_digits <- function(x) {
    if (is.null(x) || length(x) == 0 || is.na(x)) return(0L)
    sx <- sub("^-", "", as.character(x[1L]))
    if (!grepl("\\.", sx)) return(0L)
    nchar(sub("^[^.]*\\.", "", sx))
  }
  
  # Helper: apply a +/- half-up adjustment to a value
  adjust_value <- function(x, step, code) {
    if (step == 0) return(x)
    if (code == "minus") return(x - step)
    if (code == "plus")  return(x + step)
    x  # "reported"
  }
  
  # Noncentral t CI for the noncentrality parameter delta
  nct_ci <- function(t_obs, df, alpha = 0.05, max_ncp = 1000) {
    fL <- function(delta) stats::pt(t_obs, df = df, ncp = delta) - alpha / 2
    fU <- function(delta) stats::pt(t_obs, df = df, ncp = delta) - (1 - alpha / 2)
    
    lower <- -max_ncp
    upper <-  max_ncp
    
    fL_low  <- fL(lower)
    fL_high <- fL(upper)
    if (is.na(fL_low) || is.na(fL_high) || fL_low * fL_high > 0) {
      return(c(NA_real_, NA_real_))
    }
    delta_L <- uniroot(fL, lower = lower, upper = upper)$root
    
    fU_low  <- fU(upper)
    fU_high <- fU(lower)
    # safer to flip the bracket if needed
    if (is.na(fU_low) || is.na(fU_high) || fU_low * fU_high > 0) {
      fU_low  <- fU(lower)
      fU_high <- fU(upper)
      if (is.na(fU_low) || is.na(fU_high) || fU_low * fU_high > 0) {
        return(c(NA_real_, NA_real_))
      }
    }
    delta_U <- uniroot(fU, lower = lower, upper = upper)$root
    
    c(delta_L, delta_U)
  }
  
  # What information do we have?
  vals_summary <- list(m1, m2, sd1, sd2, n1, n2)
  have_summary <- all(vapply(
    vals_summary,
    function(z) !is.null(z) && length(z) > 0 && !is.na(z[1L]),
    logical(1)
  ))
  
  # Now: treat t/df as usable even without n1/n2 (for p), but they are
  # needed for converting t to d.
  have_t <- !is.null(t) && !is.null(df) &&
    !is.na(t)   && !is.na(df)
  
  if (!have_summary && !have_t) {
    stop("Provide either M/SD/N (m1, m2, sd1, sd2, n1, n2) and/or t and df.")
  }
  
  have_n_for_t <- !is.null(n1) && !is.null(n2) &&
    !is.na(n1)   && !is.na(n2)
  
  d_results <- list()
  p_results <- list()
  idx_d <- 1
  idx_p <- 1
  
  # ## 0) Append reported values rows ----
  # d_results[[idx_d]] <- data.frame(
  #   source              = "reported",
  #   direction           = NA_character_,
  #   es_type             = NA_character_,
  #   ci_method           = NA_character_,
  #   d_rounding          = NA_character_,
  #   input_adj_stats     = NA_character_,
  #   input_adj_tdf       = NA_character_,
  #   m1_used             = NA_real_,
  #   m2_used             = NA_real_,
  #   sd1_used            = NA_real_,
  #   sd2_used            = NA_real_,
  #   t_used              = NA_real_,
  #   df_used             = NA_real_,
  #   es_unrounded        = NA_real_,
  #   ci_lower_unrounded  = NA_real_,
  #   ci_upper_unrounded  = NA_real_,
  #   est_rounded         = d_est_num,
  #   ci_lower_rounded    = d_ci_lower_num,
  #   ci_upper_rounded    = d_ci_upper_num,
  #   match_est           = NA,
  #   match_ci_lower      = NA,
  #   match_ci_upper      = NA,
  #   match_all           = NA,
  #   stringsAsFactors    = FALSE
  # )
  # idx_d <- idx_d + 1
  
  # p_results[[idx_p]] <- data.frame(
  #   source              = "reported",
  #   direction           = NA_character_,
  #   p_method            = NA_character_,
  #   p_rounding          = NA_character_,
  #   input_adj_stats     = NA_character_,
  #   input_adj_tdf       = NA_character_,
  #   t_used              = NA_real_,
  #   df_used             = NA_real_,
  #   p_unrounded         = NA_real_,
  #   p_rounded           = p_est_num,
  #   match_p             = NA,
  #   stringsAsFactors    = FALSE
  # )
  # idx_p <- idx_p + 1
  
  ci_methods <- c("wald_z", "wald_t", "welch_t", "welch_z", "nct")
  p_methods  <- c("student_t", "welch_t", "student_z", "welch_z")
  
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
        
        df_s  <- n1 + n2 - 2
        N     <- n1 + n2
        sp    <- sqrt(((n1 - 1) * sd1_star^2 + (n2 - 1) * sd2_star^2) / df_s)
        if (!is.finite(sp) || sp <= 0) next
        
        # Pooled t (Student)
        t_pooled <- diff_mean / (sp * sqrt(1 / n1 + 1 / n2))
        
        # Welch variance and t
        var1 <- sd1_star^2
        var2 <- sd2_star^2
        se_welch <- sqrt(var1 / n1 + var2 / n2)
        t_welch  <- diff_mean / se_welch
        
        # Welch df
        num_w <- (var1 / n1 + var2 / n2)^2
        den_w <- (var1^2 / (n1^2 * (n1 - 1))) +
          (var2^2 / (n2^2 * (n2 - 1)))
        df_w <- num_w / den_w
        if (!is.finite(df_w) || df_w <= 0) df_w <- NA_real_
        
        # Cohen's d and Hedges' g
        d_raw <- diff_mean / sp
        J_s   <- 1 - 3 / (4 * df_s - 1)
        g_raw <- J_s * d_raw
        
        ## 1a) Effect sizes and CIs -> d_results ----
        for (es_type in c("d", "g")) {
          
          es <- if (es_type == "d") d_raw else g_raw
          
          # SE under pooled df and Welch df for Wald CIs
          se_pooled <- sqrt(N / (n1 * n2) + es^2 / (2 * df_s))
          se_welch  <- if (!is.na(df_w)) {
            sqrt(N / (n1 * n2) + es^2 / (2 * df_w))
          } else {
            NA_real_
          }
          
          for (ci_method in ci_methods) {
            
            # Compute CI for d/g
            if (ci_method == "nct") {
              # Noncentral t CI based on pooled t and df_s
              delta_ci <- nct_ci(t_pooled, df_s, alpha = alpha)
              if (any(is.na(delta_ci))) next
              
              fac_d  <- sqrt(1 / n1 + 1 / n2)
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
            
            for (d_rounding in c("half_up", "half_down")) {
              
              d_round_fun <- if (d_rounding == "half_up") {
                function(x) round_up(x, d_digits)
              } else {
                function(x) round_down(x, d_digits)
              }
              
              est_r   <- d_round_fun(es)
              lower_r <- d_round_fun(ci_lower)
              upper_r <- d_round_fun(ci_upper)
              
              d_results[[idx_d]] <- data.frame(
                source              = "summary",
                direction           = direction,
                es_type             = es_type,
                ci_method           = ci_method,
                d_rounding          = d_rounding,
                input_adj_stats     = adj_stats,
                input_adj_tdf       = NA_character_,
                m1_used             = m1_star,
                m2_used             = m2_star,
                sd1_used            = sd1_star,
                sd2_used            = sd2_star,
                t_used              = t_pooled,  # primary t used for CI (for nct)
                df_used             = df_s,
                es_unrounded        = es,
                ci_lower_unrounded  = ci_lower,
                ci_upper_unrounded  = ci_upper,
                est_rounded         = est_r,
                ci_lower_rounded    = lower_r,
                ci_upper_rounded    = upper_r,
                match_est           = if (!is.na(d_est_num)) isTRUE(all.equal(est_r,   d_est_num))      else NA,
                match_ci_lower      = if (!is.na(d_ci_lower_num)) isTRUE(all.equal(lower_r, d_ci_lower_num)) else NA,
                match_ci_upper      = if (!is.na(d_ci_upper_num)) isTRUE(all.equal(upper_r, d_ci_upper_num)) else NA,
                match_all           = if (!any(is.na(c(d_est_num, d_ci_lower_num, d_ci_upper_num)))) {
                  (est_r == d_est_num &&
                     lower_r == d_ci_lower_num &&
                     upper_r == d_ci_upper_num)
                } else NA,
                stringsAsFactors    = FALSE
              )
              
              idx_d <- idx_d + 1
            }
          }
        }
        
        ## 1b) p-values -> p_results ----
        for (p_method in p_methods) {
          
          if (p_method == "student_t") {
            t_use <- t_pooled
            df_use <- df_s
            p_unr <- 2 * (1 - stats::pt(abs(t_use), df = df_use))
          } else if (p_method == "welch_t") {
            if (is.na(df_w)) next
            t_use <- t_welch
            df_use <- df_w
            p_unr <- 2 * (1 - stats::pt(abs(t_use), df = df_use))
          } else if (p_method == "student_z") {
            t_use <- t_pooled
            df_use <- Inf
            p_unr <- 2 * (1 - stats::pnorm(abs(t_use)))
          } else { # welch_z
            t_use <- t_welch
            df_use <- Inf
            p_unr <- 2 * (1 - stats::pnorm(abs(t_use)))
          }
          
          for (p_rounding in c("half_up", "half_down")) {
            
            p_round_fun <- if (p_rounding == "half_up") {
              function(x) round_up(x, p_digits)
            } else {
              function(x) round_down(x, p_digits)
            }
            
            p_rounded <- p_round_fun(p_unr)
            
            p_results[[idx_p]] <- data.frame(
              source              = "summary",
              direction           = direction,
              p_method            = p_method,
              p_rounding          = p_rounding,
              input_adj_stats     = adj_stats,
              input_adj_tdf       = NA_character_,
              t_used              = t_use,
              df_used             = df_use,
              p_unrounded         = p_unr,
              p_rounded           = p_rounded,
              match_p             = if (!is.na(p_est_num)) isTRUE(all.equal(p_rounded, p_est_num)) else NA,
              stringsAsFactors    = FALSE
            )
            
            idx_p <- idx_p + 1
          }
        }
      }
    }
  }
  
  #### 2) From t statistic and df ####
  if (have_t) {
    t_val <- as.numeric(t)
    df_t  <- as.numeric(df)
    
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
      
      ## 2a) Effect sizes from t: only if n1 and n2 available
      if (have_n_for_t) {
        N     <- n1 + n2
        df_s  <- n1 + n2 - 2
        fac_d <- sqrt(1 / n1 + 1 / n2)
        d_from_t <- t_star * fac_d
        
        # Allow possible sign flip (authors may have reversed group labels)
        for (sign_flip in c(1, -1)) {
          
          es_base  <- d_from_t * sign_flip
          direction <- if (sign_flip == 1) "t_sign" else "neg_t_sign"
          
          # Hedges' J for this df_star
          J_t   <- 1 - 3 / (4 * df_star - 1)
          g_base <- J_t * es_base
          
          for (es_type in c("d", "g")) {
            
            es <- if (es_type == "d") es_base else g_base
            
            # SE under pooled df_s and df_star
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
              
              for (d_rounding in c("half_up", "half_down")) {
                
                d_round_fun <- if (d_rounding == "half_up") {
                  function(x) round_up(x, d_digits)
                } else {
                  function(x) round_down(x, d_digits)
                }
                
                est_r   <- d_round_fun(es)
                lower_r <- d_round_fun(ci_lower)
                upper_r <- d_round_fun(ci_upper)
                
                d_results[[idx_d]] <- data.frame(
                  source              = "t_df",
                  direction           = direction,
                  es_type             = es_type,
                  ci_method           = ci_method,
                  d_rounding          = d_rounding,
                  input_adj_stats     = NA_character_,
                  input_adj_tdf       = adj_tdf,
                  m1_used             = NA_real_,
                  m2_used             = NA_real_,
                  sd1_used            = NA_real_,
                  sd2_used            = NA_real_,
                  t_used              = t_star,
                  df_used             = df_star,
                  es_unrounded        = es,
                  ci_lower_unrounded  = ci_lower,
                  ci_upper_unrounded  = ci_upper,
                  est_rounded         = est_r,
                  ci_lower_rounded    = lower_r,
                  ci_upper_rounded    = upper_r,
                  match_est           = if (!is.na(d_est_num)) isTRUE(all.equal(est_r,   d_est_num))      else NA,
                  match_ci_lower      = if (!is.na(d_ci_lower_num)) isTRUE(all.equal(lower_r, d_ci_lower_num)) else NA,
                  match_ci_upper      = if (!is.na(d_ci_upper_num)) isTRUE(all.equal(upper_r, d_ci_upper_num)) else NA,
                  match_all           = if (!any(is.na(c(d_est_num, d_ci_lower_num, d_ci_upper_num)))) {
                    (est_r == d_est_num &&
                       lower_r == d_ci_lower_num &&
                       upper_r == d_ci_upper_num)
                  } else NA,
                  stringsAsFactors    = FALSE
                )
                
                idx_d <- idx_d + 1
              }
            }
          }
        }
      } # end have_n_for_t
      
      ## 2b) p-values -> p_results (works even without n1/n2) ----
      for (p_method in p_methods) {
        
        if (p_method %in% c("student_t", "welch_t")) {
          # Without n1/n2 we cannot distinguish, use df_star in both cases.
          t_use <- t_star
          df_use <- df_star
          p_unr <- 2 * (1 - stats::pt(abs(t_use), df = df_use))
        } else { # student_z, welch_z: treat t_star as z
          t_use <- t_star
          df_use <- Inf
          p_unr <- 2 * (1 - stats::pnorm(abs(t_use)))
        }
        
        for (p_rounding in c("half_up", "half_down")) {
          
          p_round_fun <- if (p_rounding == "half_up") {
            function(x) round_up(x, p_digits)
          } else {
            function(x) round_down(x, p_digits)
          }
          
          p_rounded <- p_round_fun(p_unr)
          
          p_results[[idx_p]] <- data.frame(
            source              = "t_df",
            direction           = if (have_n_for_t) "t_sign" else NA_character_,
            p_method            = p_method,
            p_rounding          = p_rounding,
            input_adj_stats     = NA_character_,
            input_adj_tdf       = adj_tdf,
            t_used              = t_use,
            df_used             = df_use,
            p_unrounded         = p_unr,
            p_rounded           = p_rounded,
            match_p             = if (!is.na(p_est_num)) isTRUE(all.equal(p_rounded, p_est_num)) else NA,
            stringsAsFactors    = FALSE
          )
          
          idx_p <- idx_p + 1
        }
      }
    }
  }
  
  d_out <- do.call(rbind, d_results)
  p_out <- do.call(rbind, p_results)
  
  # Order d_out by how well it matches d and CIs (reported row unaffected: match_all = NA)
  d_out <- d_out[order(!d_out$match_all,
                       !d_out$match_est,
                       !d_out$match_ci_lower,
                       !d_out$match_ci_upper), ]
  
  list(
    reproduced = any(d_out$match_all, na.rm = TRUE),
    d_results  = d_out,
    p_results  = p_out
  )
}

# example
res <- identify_d_method(
  n1 = 50, 
  n2 = 48,
  m1 = 10.3,
  m2 = 8.7,
  sd1 = 3.1,
  sd2 = 2.8,
  t = 1.34,
  df = 97,
  d_est = 0.55,
  d_ci_lower = 0.20,
  d_ci_upper = 0.90,
  d_digits = 2
)

res_d <- res$d_results |>
  as_tibble() |>
  mutate(original = source=="reported") |>
  arrange(est_rounded, original) |>
  rownames_to_column(var = "rowname") |>
  mutate(rowname = as.numeric(as.character(rowname)))
  
res_d |>
  filter(source != "reported" & direction %in% c("m1_minus_m2", "t_sign")) |>
  ggplot(aes(y = rowname, 
             x = est_rounded, xmin = ci_lower_rounded, xmax = ci_upper_rounded,
             color = original)) +
  ggstance::geom_linerangeh() +
  geom_point() +
  facet_wrap(~ source)

res_d |>
  filter(source != "reported" & direction %in% c("m1_minus_m2", "t_sign")) |>
  summarize(min_est_rounded = min(est_rounded),
            max_est_rounded = max(est_rounded))


res_p <- res$p_results |>
  as_tibble() |>
  mutate(original = source=="reported") |>
  arrange(p_unrounded, original) |>
  rownames_to_column(var = "rowname") |>
  mutate(rowname = as.numeric(as.character(rowname)))

res_p |>
  filter(source != "reported" & direction %in% c("m1_minus_m2", "t_sign")) |>
  ggplot(aes(y = rowname, 
             x = p_unrounded,
             color = original)) +
  geom_point()

res_p |>
  filter(source != "reported" & direction %in% c("m1_minus_m2", "t_sign")) |>
  summarize(min_p_unrounded = min(p_unrounded),
            max_p_unrounded = max(p_unrounded))



  

