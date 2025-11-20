multiverse_one_or_two_way_between_anova <- function(
    m,          # means: vector (1-way) or matrix (2-way)
    sd,         # sds:   same shape as m
    n,          # ns:    same shape as m, or scalar/vec recycled
    # reported values to match (for a given effect)
    F_est      = NULL,
    df1_est    = NULL,
    df2_est    = NULL,
    p_est      = NULL,
    peta2_est  = NULL,
    F_digits   = 2,
    p_digits   = 3,
    peta2_digits = 3,
    alpha_peta2 = 0.10,         # CI level for partial eta2
    
    # researcher degrees of freedom
    unbiased_options = c(TRUE, FALSE),   # biased vs unbiased SD
    input_adj_codes  = c("reported", "minus", "plus"),  # rounding of M/SD
    F_rounding_set   = c("half_up", "half_down", "bankers"),
    p_rounding_set   = c("half_up", "half_down", "bankers"),
    peta2_rounding_set = c("half_up", "half_down", "bankers")
) {
  require(roundwork)
  require(dplyr)
  
  ## ---- validation / coerce reported values ----
  F_est_num     <- if (!is.null(F_est))     as.numeric(F_est)     else NA_real_
  df1_est_num   <- if (!is.null(df1_est))   as.numeric(df1_est)   else NA_real_
  df2_est_num   <- if (!is.null(df2_est))   as.numeric(df2_est)   else NA_real_
  p_est_num     <- if (!is.null(p_est))     as.numeric(p_est)     else NA_real_
  peta2_est_num <- if (!is.null(peta2_est)) as.numeric(peta2_est) else NA_real_
  
  # allowed sets
  allowed_F_round <- c("half_up", "half_down", "bankers")
  allowed_p_round <- c("half_up", "half_down", "bankers")
  allowed_eta_round <- c("half_up", "half_down", "bankers")
  allowed_adj <- c("reported", "minus", "plus")
  
  if (length(bad <- setdiff(F_rounding_set, allowed_F_round))) {
    stop("F_rounding_set contains invalid values: ",
         paste0(bad, collapse = ", "),
         ". Allowed: ", paste0(allowed_F_round, collapse = ", "))
  }
  if (length(bad <- setdiff(p_rounding_set, allowed_p_round))) {
    stop("p_rounding_set contains invalid values: ",
         paste0(bad, collapse = ", "),
         ". Allowed: ", paste0(allowed_p_round, collapse = ", "))
  }
  if (length(bad <- setdiff(peta2_rounding_set, allowed_eta_round))) {
    stop("peta2_rounding_set contains invalid values: ",
         paste0(bad, collapse = ", "),
         ". Allowed: ", paste0(allowed_eta_round, collapse = ", "))
  }
  if (length(bad <- setdiff(input_adj_codes, allowed_adj))) {
    stop("input_adj_codes contains invalid values: ",
         paste0(bad, collapse = ", "),
         ". Allowed: ", paste0(allowed_adj, collapse = ", "))
  }
  
  # helper: digits in a number
  get_digits <- function(x) {
    if (is.null(x) || length(x) == 0 || is.na(x)) return(0L)
    sx <- sub("^-", "", as.character(x[1L]))
    if (!grepl("\\.", sx)) return(0L)
    nchar(sub("^[^.]*\\.", "", sx))
  }
  
  # helper: rounding functions
  make_round_fun <- function(which, digits) {
    switch(
      which,
      "half_up"   = function(x) round_up(x,   digits),
      "half_down" = function(x) round_down(x, digits),
      "bankers"   = function(x) round(x,      digits),
      stop("Unknown rounding mode: ", which)
    )
  }
  
  # coerce m, sd, n to matrices
  m_mat  <- if (is.matrix(m)) m else matrix(m, ncol = 1)
  sd_mat <- if (is.matrix(sd)) sd else matrix(sd, ncol = ncol(m_mat))
  if (is.vector(n)) {
    # recycle n into matrix
    n_mat <- matrix(n, nrow = nrow(m_mat), ncol = ncol(m_mat))
  } else {
    n_mat <- if (is.matrix(n)) n else matrix(n, ncol = ncol(m_mat))
  }
  
  # global decimal assumption for means/sds (like your ANOVA helper)
  dig_m  <- max(vapply(as.numeric(m_mat),  get_digits, integer(1)))
  dig_sd <- max(vapply(as.numeric(sd_mat), get_digits, integer(1)))
  step_m  <- 0.5 * 10^(-dig_m)
  step_sd <- 0.5 * 10^(-dig_sd)
  
  # helper: adjust all cells together
  adjust_matrix <- function(x, step, code) {
    if (code == "minus") return(x - step)
    if (code == "plus")  return(x + step)
    x  # "reported"
  }
  
  results <- list()
  idx <- 1
  
  for (unbiased in unbiased_options) {
    for (adj_code in input_adj_codes) {
      
      m_star  <- adjust_matrix(m_mat,  step_m,  adj_code)
      sd_star <- adjust_matrix(sd_mat, step_sd, adj_code)
      
      # ind_twoway_second returns the full ANOVA table and partial eta2/CI
      fit <- ind_twoway_second(
        m        = m_star,
        sd       = sd_star,
        n        = n_mat,
        unbiased = unbiased,
        sig_level = alpha_peta2
      )
      
      anova_tab <- fit$anova.table
      # partial eta2 / CI for row, col, interaction
      eta_tab   <- fit$omnibus.es
      
      # Extract within df (error df)
      df_within <- anova_tab["Within", "df"]
      
      # We care about each between effect separately
      effect_rows <- c("Between (row)", "Between (col)", "Between (row * col)")
      effect_rows <- effect_rows[effect_rows %in% rownames(anova_tab)]
      
      for (effect in effect_rows) {
        F_val  <- anova_tab[effect, "F"]
        df1    <- anova_tab[effect, "df"]
        df2    <- df_within
        
        # p and partial eta2 for this effect
        p_unr     <- 1 - stats::pf(F_val, df1, df2)
        peta2     <- eta_tab[effect, "partial_etasq"]
        peta2_L   <- eta_tab[effect, "partial_etasq_lower"]
        peta2_U   <- eta_tab[effect, "partial_etasq_upper"]
        
        for (F_round in F_rounding_set) {
          F_round_fun <- make_round_fun(F_round, F_digits)
          F_rep <- F_round_fun(F_val)
          
          for (p_round in p_rounding_set) {
            p_round_fun <- make_round_fun(p_round, p_digits)
            p_rep <- p_round_fun(p_unr)
            
            for (eta_round in peta2_rounding_set) {
              eta_round_fun <- make_round_fun(eta_round, peta2_digits)
              eta_rep   <- eta_round_fun(peta2)
              etaL_rep  <- eta_round_fun(peta2_L)
              etaU_rep  <- eta_round_fun(peta2_U)
              
              results[[idx]] <- data.frame(
                effect          = effect,
                unbiased        = unbiased,
                input_adj_code  = adj_code,
                F_rounding      = F_round,
                p_rounding      = p_round,
                peta2_rounding  = eta_round,
                
                F_unrounded     = F_val,
                df1_used        = df1,
                df2_used        = df2,
                p_unrounded     = p_unr,
                peta2_unrounded = peta2,
                peta2_ci_lower_unrounded = peta2_L,
                peta2_ci_upper_unrounded = peta2_U,
                
                F_rounded       = F_rep,
                p_rounded       = p_rep,
                peta2_rounded   = eta_rep,
                peta2_ci_lower_rounded = etaL_rep,
                peta2_ci_upper_rounded = etaU_rep,
                
                match_F     = if (!is.na(F_est_num))     isTRUE(all.equal(F_rep,     F_est_num))     else NA,
                match_df1   = if (!is.na(df1_est_num))   isTRUE(all.equal(df1,       df1_est_num))   else NA,
                match_df2   = if (!is.na(df2_est_num))   isTRUE(all.equal(df2,       df2_est_num))   else NA,
                match_p     = if (!is.na(p_est_num))     isTRUE(all.equal(p_rep,     p_est_num))     else NA,
                match_peta2 = if (!is.na(peta2_est_num)) isTRUE(all.equal(eta_rep,   peta2_est_num)) else NA,
                stringsAsFactors = FALSE
              )
              
              idx <- idx + 1
            }
          }
        }
      }
    }
  }
  
  out <- dplyr::bind_rows(results)
  
  # Order by best match (if reported values supplied)
  if (!all(is.na(out$match_F)) ||
      !all(is.na(out$match_p)) ||
      !all(is.na(out$match_peta2))) {
    out <- out[order(
      !out$match_F,
      !out$match_df1,
      !out$match_df2,
      !out$match_p,
      !out$match_peta2
    ), ]
  }
  
  out
}


labels <- c("groups 1 (A vs B)", "groups 2 (C vs D)", "groups 1 x groups 2 interaction")

m <- matrix(c(5.00, 2.69, 
              4.83, 5.54), 
            ncol = 2)

sd <- matrix(c(2.99, 2.57, 
               2.71, 1.84), 
             ncol = 2)

n <- matrix(c(40, 20, 
              35, 10), 
            ncol = 2)

res_anova <- multiverse_one_or_two_way_between_anova(
  m = m, 
  sd = sd, 
  n = n,
  F_est     = 4.32,   # if you have a reported F, otherwise leave NULL
  df1_est   = 1,
  df2_est   = 99,
  p_est     = 0.041,
  peta2_est = 0.065,
  F_digits  = 2,
  p_digits  = 3,
  peta2_digits = 2
)

head(res_anova)

res_sum <- res_anova |>
  group_by(effect) |>
  summarize(min_p_rounded = min(p_rounded),
            max_p_rounded = max(p_rounded),
            min_peta2_rounded = min(peta2_rounded),
            max_peta2_rounded = max(peta2_rounded))
