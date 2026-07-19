library(shiny)
library(rhandsontable)

# main function -----------------------------------------------------------

base_p_cat_format <- function(df) {
  
  a <- df
  
  a <- a %>% dplyr::rename(study = Study, var= Variable) 
  
  colnames (a) [tolower(colnames(a)) == "variable"] <- "var"
  colnames (a) [startsWith(colnames (a), "N")] <- paste0("t", 
                                                         substr(colnames (a) [startsWith(colnames (a), "N")],2,nchar(colnames(a))))
  colnames (a) <- lapply(colnames (a), tolower)
  
  if ("stat" %in% colnames (a)) {a <- a %>% dplyr::group_by(study, var) %>% tidyr::fill(stat, .direction = "updown")}
  
  if ("p" %in% colnames (a)) {a <- a %>% dplyr::group_by(study, var) %>% tidyr::fill(p, .direction = "updown")}
  
  a <- a %>% dplyr::group_by(study, var) %>% tidyr::fill(p, starts_with("t"), .direction = "updown") %>% dplyr::ungroup()
  
  #keep stat till the end
  stat = 0
  
  if ("stat" %in% colnames(a)) {
    if (any(a %>% dplyr::group_by(study, var) %>% 
            dplyr::summarise(n = dplyr::n_distinct(stat), .groups = "drop") %>% 
            dplyr::pull(n) > 1)) {
      return("Data error- some 'stat' variables don't match for different levels")
    }
  }
  
  if ("stat" %in% colnames(a)) {
    stat_df <- a %>% dplyr::distinct(study, var, stat)
  } else {stat <- 1}
  
  a$stat <- NULL
  
  a <- a %>% dplyr::mutate(across(starts_with(c("n", "t")), as.numeric))  
  
  #remove and add temporary separator
  names <- colnames(a) [substr(colnames(a),1,1) %in% c("n","t")]
  colnames (a) [colnames (a) %in% names] <- paste0(substr(names,1,1),"_", gsub("\\D", "", names))
  
  #add p if none
  if (!"p" %in% colnames (a)) {a$p <- NA}
  
  #get names;
  nm1 <- colnames(a)[startsWith(colnames(a), "n")]
  nm2 <- colnames(a)[startsWith(colnames(a), "t")]
  
  #check groups
  if (!"group" %in% colnames(a)) {a$group <- rowSums(!is.na( a [,nm1]))}
  
  #stop if not matching#
  if (length(nm1) != length(nm2) | any(a$group != rowSums(!is.na(a[, nm2])))) {
    return ("Data error- some data for groups are missing or group numbers don't match")
  }
  
  #stop if blank rows
  if (any(rowSums(is.na(a[, c(nm1, nm2)])) == length(c(nm1,nm2)))) {
    return ("Data error- some rows are blank")}
  
  #number of levels
  # Create or fix level column
  
  if ("level" %in% colnames(a)) {a$level <- as.numeric(a$level)}
  if ("level" %in% colnames(a) && any(is.na(a$level))) {
    # If any level is NA for a study-var, set all levels to NA for that study-var
    a <- a %>%
      dplyr::group_by(study, var) %>%
      dplyr::mutate(level = if(any(is.na(level))) NA else level) %>%
      dplyr::ungroup()
  }
  
  if (!"level" %in% colnames(a)) {a$level <- NA}
  
  # Assign levels based on total n (descending) where missing
  if (!"calculated_level" %in% colnames(a)) {a$calculated_level <- NA_character_}
  a <- a %>% dplyr::mutate(
    calculated_level = dplyr::case_when(!is.na(calculated_level) & 
                                   startsWith(tolower(calculated_level), "y") ~ "yes",
             .default = NA_character_))

  a <- a %>%
    dplyr::mutate(total_n = rowSums(dplyr::across(all_of(nm1)), na.rm = TRUE)) %>%
    dplyr::arrange(study, var, calculated_level, desc(total_n)) %>% dplyr::group_by(study, var) %>%
    dplyr::mutate(level = ifelse(is.na(level), dplyr::row_number(), level)) %>%
    dplyr::select(-total_n) %>% dplyr::mutate(n_all = dplyr::n()) %>% dplyr::ungroup() 
  
  # Create calculated level
  if (any(a$n_all == 1)) {
    
    l1 <- a %>% dplyr::filter(n_all == 1) %>% dplyr::mutate(n_all = 2, level=2) #level 2 calc level =1
    l2 <- l1 %>% dplyr::mutate(level = 1, calculated_level = "yes")
    
    for (i in seq_along(nm1)) {l2[[nm1[i]]] <- l2[[nm2[i]]] - l2[[nm1[i]]]} #subtract l2 from l1
    
    a <- dplyr::bind_rows(a %>% dplyr::filter(n_all > 1), l1, l2)
  }
  
  if (a %>% dplyr::group_by(study, var) %>% 
      dplyr::summarise(n = dplyr::n_distinct(calculated_level, na.rm = TRUE), .groups = "drop") %>%
      dplyr::filter(n > 1) %>% nrow() > 0) {
    return("Some study-var combinations have multiple calculated_level values")}
  
  
  #check data matches
  # Calculate correct sums grouped by study and var
  chk <- a %>% dplyr::group_by(study, var) %>%
    dplyr::summarise(across(all_of(nm1), ~ sum(., na.rm = TRUE)), .groups = "drop")
  
  a <- a %>%
    dplyr::left_join(chk, by = c("study", "var"), suffix = c("", "_correct")) %>%
    dplyr::mutate(N_flag = {
      m <- character(dplyr::n())
      for (i in seq_along(nm2)) {
        n_col <- paste0(sub("t", "n", nm2[i]), "_correct")
        m1 <- get(nm2[i]) != get(n_col)
        m[m1] <- paste0(m[m1], ifelse(m[m1] == "", "", "; "), sub("t_", "n", nm2[i]), " ≠ ", sub("t_", "N", nm2[i]))}
      m},
      dplyr::across(all_of(nm2), ~ get(paste0(sub("t", "n", dplyr::cur_column()), "_correct")))) %>%
    dplyr::select(-ends_with("_correct"))
  
  #now arrange into standard format
  if (all(is.na(a$p))) {a$p <- NULL}   #remove p-value column if null
  
  nm <- c("study","var","level","n_all","group", nm1, nm2)
  if ("p" %in% colnames(a)) {nm <- c(nm, "p")}
  if ("stat" %in% colnames(a)) {nm <- c(nm, "stat")}
  nm <- c(nm, "N_flag", "calculated_level")
  a <- a [,nm]
  
  colnames(a)[4] <- "levels_no"
  colnames(a)[startsWith(colnames(a), "n")] <- paste0("n", gsub("\\D", "", colnames(a)[startsWith(colnames(a), "n")]))
  colnames(a)[startsWith(colnames(a), "t")] <- paste0("N", gsub("\\D", "", colnames(a)[startsWith(colnames(a), "t")]))
  
  #check p-vals/tests are correct
  if ("p" %in% colnames(a)) {
    p_check <- a %>% dplyr::group_by(study, var) %>%
      dplyr::summarise(n_unique_p = dplyr::n_distinct(p), .groups = "drop") %>%
      dplyr::filter(n_unique_p > 1)
    
    if (nrow(p_check) > 0) {
      return("Data error- some p-values don't match for different levels")
    }
  }
  
  # Add p-value to each level
  if ("p" %in% colnames(a)) {
    a <- a %>% dplyr::group_by(study, var) %>% dplyr::mutate(p = dplyr::first(p[!is.na(p)])) %>%
      dplyr::ungroup() }
  
  #add in stat
  if (stat == 0) {a <- dplyr::left_join(a, stat_df, by = c("study", "var"))}
  
  a <- a %>% dplyr::arrange (as.numeric(gsub("\\D", "", study)), var, level)
  
  if (any(is.na(c(a$n1, a$n2, a$N1, a$N2)))) {
    return("Error- some data are missing for numbers in n1/n2/N1/N2") }
  
  return (a)
}

base_p_cat_calcs <- function(a) {
  b <- a %>%
    dplyr::select(study, var, level, starts_with("n"), -N_flag) %>%
    tidyr::pivot_longer(cols = starts_with("n"), names_to = "k", values_to = "v") %>%
    
    dplyr::mutate(pre = gsub("[0-9]+", "", k), num = gsub("\\D", "", k)) %>%  # Split column names into prefix and group number
    dplyr::arrange(level, num) %>%
    tidyr::unite("k", pre, level, num, sep = "_") %>%
    dplyr::mutate(k = factor(k, levels = unique(k))) %>% 
    tidyr::pivot_wider(names_from = k, values_from = v)
  
  b <- dplyr::left_join(b, a %>% dplyr::select(study, var, levels_no ,group) %>%
                          dplyr::distinct(study, var, .keep_all = TRUE), by = c("study", "var")) %>% 
    dplyr::select(-starts_with("N", ignore.case = FALSE))
  
  
  p_cols <- c("p_chi", "p_chi_yates", "p_fisher", "p_lr", "p_cmh", "p_midp_exact", "p_midp_classic")
  b[p_cols] <- NA
  b$chi_warning <- NA
  n_cols <- which(startsWith(colnames(b), "n"))
  
  #functions for tests
  #chisq
  chi <- function(x, correct = FALSE) {
    
    if (nrow(x) == 2L && ncol(x) == 2L) { #fast 2*2
      a <- x[1,1]; b <- x[1,2]; c <- x[2,1]; d <- x[2,2]; n <- a + b + c + d
      
      if (correct) {chi_stat <- n * (max(abs(a*d - b*c) - n/2, 0))^2 / ((a+b)*(c+d)*(a+c)*(b+d))
      } else {chi_stat <- n * (a*d - b*c)^2 / ((a+b)*(c+d)*(a+c)*(b+d))}
      
      p <- pchisq(chi_stat, 1L, lower.tail = FALSE)
      
      row1 <- a + b; row2 <- c + d; col1 <- a + c; col2 <- b + d
      exp_min <- min(row1*col1, row1*col2, row2*col1, row2*col2) / n
      warning_msg <- if (exp_min < 5) "Expected cells warning" else ""
      
      return(list(p = p, warning = warning_msg))
      
    } else {  #fast >2*2
      nr <- nrow(x); nc <- ncol(x)
      row <- .rowSums(x, nr, nc); col <- .colSums(x, nr, nc); n <- sum(row)
      exp <- tcrossprod(row, col) / n
      
      if (correct && nr == 2L && nc == 2L) {chi_stat <- sum((pmax(abs(x - exp) - 0.5, 0))^2 / exp)
      } else {diff <- x - exp; chi_stat <- sum(diff * diff / exp)}
      
      df <- (nr - 1L) * (nc - 1L)
      p <- pchisq(chi_stat, df, lower.tail = FALSE)
      
      warning_msg <- if (length(exp) == 4L && any(exp < 5) || 
                         length(exp) > 4L && (any(exp < 1) || sum(exp < 5) > 0.2 * length(exp))) { "Expected cells warning"
      } else ""
      
      return(list(p = p, warning = warning_msg))
    }
  }
  
  fisher_fast <- function(x, alternative = "two.sided", B = 10000) {
    nr <- nrow(x)
    nc <- ncol(x)
    
    # 2x2 tables - exact calculation
    if (nr == 2L && nc == 2L) {
      a <- x[1,1]; b <- x[1,2]; c <- x[2,1]; d <- x[2,2]
      
      m <- a + c; n <- b + d; k <- a + b
      
      if (alternative == "two.sided") { # Two-sided p-value
        p_obs <- dhyper(a, m, n, k)
        lo <- max(0L, k - n); hi <- min(k, m); support <- lo:hi
        p_all <- dhyper(support, m, n, k)
        
        # Sum probabilities <= observed
        p_value <- sum(p_all[p_all <= p_obs * (1 + 1e-7)])
      } else if (alternative == "less") {
        p_value <- phyper(a, m, n, k)
      } else if (alternative == "greater") {
        p_value <- phyper(a - 1, m, n, k, lower.tail = FALSE)
      } 
      
      return(p_value)
    }
    
    # Larger tables - simulation only (alternative ignored)
    sr <- rowSums(x); sc <- colSums(x)
    x <- x[sr > 0, sc > 0, drop = FALSE]
    
    # Calculate test statistic
    STATISTIC <- -sum(lfactorial(x))
    
    # Simulate using C function
    tmp <- .Call(stats:::C_Fisher_sim, rowSums(x), colSums(x), as.integer(B))
    
    # Calculate p-value
    almost.1 <- 1 + 64 * .Machine$double.eps
    p_value <- (1 + sum(tmp <= STATISTIC/almost.1)) / (B + 1)
    
    return(max(0, min(1, p_value)))
  }    
  
  # Likelihood ratio test
  lr <- function(obs) {
    row <- rowSums(obs)
    col <- colSums(obs)
    n <- sum(obs)
    exp <- outer(row, col) / n
    
    #calculate G-stat (lr chi-square), G² = 2 * sum(obs * log(obs/exp)), Only cells obs> 0 to avoid log(0)
    g_stat <- 2 * sum(obs[obs > 0] * log(obs[obs > 0] / exp[obs > 0]))
    df <- (nrow(obs) - 1) * (ncol(obs) - 1)
    p <- pchisq(g_stat, df, lower.tail = FALSE)
    return(p)
  }
  
  # Calculate CMH general association test for a single 2D table
  cmh <- function(x) {
    R <- nrow(x); C <- ncol(x)
    
    # Total and marginal probabilities
    nt <- sum(x)
    pr <- rowSums(x) / nt
    pc <- colSums(x) / nt
    
    # Expected frequencies (as vector)
    m <- as.vector(nt * outer(pr, pc))
    n <- as.vector(x)
    
    # Covariance matrix
    V1 <- diag(pr) - pr %*% t(pr)
    V2 <- diag(pc) - pc %*% t(pc)
    V <- (nt^2 / (nt - 1)) * kronecker(V2, V1)
    
    # Contrast matrix for general association
    # A = kronecker of reduced identity matrices
    A_row <- cbind(diag(R - 1), rep(0, R - 1))
    A_col <- cbind(diag(C - 1), rep(0, C - 1))
    A <- kronecker(A_col, A_row)
    
    # CMH test statistic
    AVA <- A %*% V %*% t(A)
    Q <- as.numeric(t(n - m) %*% t(A) %*% solve(AVA) %*% A %*% (n - m))
    
    # Degrees of freedom for general association
    df <- (R - 1) * (C - 1)
    
    # P-value
    p_value <- pchisq(Q, df, lower.tail = FALSE)
    
    return(p_value)
  }
  
  # Calculate mid-p exact test for 2x2 table
  midp_exact <- function(x) {
    a1 <- x[1, 1]  # exposed, outcome
    b1 <- x[1, 2]  # exposed, no outcome
    a0 <- x[2, 1]  # unexposed, outcome
    b0 <- x[2, 2]  # unexposed, no outcome
    
    # Fisher's exact test (one-sided)
    lteq <- fisher_fast(x, alternative = "less")
    gteq <- fisher_fast(x, alternative = "greater")
    
    # Mid-p calculation
    pval1 <- 0.5 * (lteq - gteq + 1)
    one_sided <- min(pval1, 1 - pval1)
    two_sided <- 2 * one_sided
    
    return(two_sided)
  }
  
  # Helper function: Calculate SAS-style mid-p value for 2x2 table
  midp_classic <- function(x, p_fisher) {
    # Extract cell counts
    a <- x[1, 1]; b <- x[1, 2]
    c <- x[2, 1]; d <- x[2, 2]
    n <- sum(x)
    
    # Calculate hypergeometric probability for observed table
    # P = (a+b)! * (c+d)! * (a+c)! * (b+d)! / (n! * a! * b! * c! * d!)
    
    # Use log-gamma for numerical stability
    log_prob <- (lgamma(a + b + 1) + lgamma(c + d + 1) + 
                   lgamma(a + c + 1) + lgamma(b + d + 1) -
                   lgamma(n + 1) - lgamma(a + 1) - lgamma(b + 1) - 
                   lgamma(c + 1) - lgamma(d + 1))
    
    prob <- exp(log_prob)
    
    # Mid-p = Fisher's exact p-value - 0.5 * P(observed)
    midp <- p_fisher - 0.5 * prob
    
    return(midp)
  }
  
  for (i in seq_len(nrow(b))) {
    #make matrix of values
    counts <- na.omit(as.numeric(b[i, n_cols]))
    x <- matrix(counts, byrow = TRUE, nrow = b$levels_no[i])
    
    # Chi-square test (uncorrected)
    chisq <- chi(x, correct = FALSE)
    b$p_chi[i] <- chisq$p
    b$chi_warning[i] <- chisq$warning
    
    b$p_chi_yates[i] <- chi(x, correct = TRUE)$p 
    
    b$p_fisher[i] <- fisher_fast(x, B = 10000)
    
    b$p_lr[i] <- lr(x) #likelihood ratio
    
    b$p_cmh[i] <- cmh(x) #cmh test
    
    if (b$group[i] == 2 && b$levels_no[i] == 2) { 
      b$p_midp_exact[i] <- midp_exact(x) # mid-p- see epitools, Rodman 1988
      b$p_midp_classic[i] <- midp_classic(x, b$p_fisher[i])
    }
    
  }
  
  #now merge original p/test back in
  if (!"stat" %in% colnames(a)) {a$stat <- NA}
  a$stat <- tolower(a$stat)
  a$stat1 <- a$stat #keep originals
  if (!"p" %in% colnames(a)) {a$p <- NA}
  
  b <- dplyr::left_join(b, a %>% dplyr::select(study, var, p, stat, stat1) %>% 
                          dplyr::distinct(study, var, .keep_all = TRUE), by = c("study", "var"))
  
  #get pvalues into usuable form
  b <- b %>% dplyr::mutate(p_n = suppressWarnings(as.numeric(p)),
                           p_g = dplyr::case_when(grepl("ns", p, ignore.case = TRUE) ~ 0.05,
                                                  grepl(">", p) ~ suppressWarnings(as.numeric(sub('.*>\\s*', '', p))), .default = NA),
                           p_l = dplyr::case_when(grepl("<", p) ~ suppressWarnings(as.numeric(sub('.*<\\s*', '', p))),
                                                  .default = NA), p = ifelse(is.na(p_n), tolower(p), as.character(p_n)))
  
  #add digits
  b$p_d <- ifelse(!is.na(b$p_n),
                  ifelse(abs(b$p_n - round(b$p_n)) > .Machine$double.eps^0.5,
                         nchar(sub('^\\d+\\.', '', sub('0+$', '', as.character(abs(b$p_n))))), 0), NA)
  
  #round all p-values to # digits
  b [, p_cols] <- lapply(b[, p_cols], function(x) {
    ifelse(!is.na(b$p_n), round(x, b$p_d), round(x, pmax(b$p_d, 3, na.rm = TRUE)))})
  
  max_digits = ifelse(sum(!is.na(b$p_d))>0, max(b$p_d, na.rm = T), 3)

  #compare p results
  nm_p <- c("p_chi", "p_fisher", "p_chi_yates", "p_lr", "p_cmh", "p_midp_exact", "p_midp_classic")
  names (nm_p) <- c("chisq", "fisher", "chisqc", "lr", "mh", "midp","midp")
  
  b$p_match1 <- sapply(seq_len(nrow(b)), function(i) {
    col <- nm_p[b$stat[i]]
    if (is.na(col)) return(NA)
    b[[col]] [i]
  })
  b$p_match2 <- ifelse(!is.na(b$stat) & b$stat == "midp", b$p_midp_classic, NA)
  
  #check against thresholds
  # check if any calculated p-value meets threshold
  thresh <- function(p, thrsh, comp = ">") {
    matches <- if (comp == ">") p > thrsh else p < thrsh
    any(matches, na.rm = TRUE)}
  
  # best match test
  test <- function(p, thrsh, comp = ">") {
    matches <- if (comp == ">") p > thrsh else p < thrsh
    if (!any(matches, na.rm = TRUE)) return(NA)
    names(nm_p)[which.max(matches)]}
  
  b$m_t_comment <- b$m_t <-  NA
  for (i in seq_len(nrow(b))) {
    if (is.na(b$p_g[i]) && is.na(b$p_l[i])) next # no thresholds so skip
    p_vals <- as.numeric(b[i, nm_p])
    table <- ifelse(b$group[i] * b$levels_no[i] <= 4, "2x2", ">2x2")
    
    #No reported test
    if (is.na(b$stat[i])) {
      if (!is.na(b$p_g[i])) {# Determine which threshold to use (prefer p_g if both exist)
        thrsh <- b$p_g[i]
        comp <- ">"
      } else {
        thrsh <- b$p_l[i]
        comp <- "<"
      }
      
      if (thresh(p_vals, thrsh, comp)) {
        best <- test(p_vals, thrsh, comp)
        b$m_t[i] <- "Match"
        b$m_t_comment[i] <- paste0("No reported test, ", table, ", test=", best)
      } else {
        b$m_t[i] <- "No match"
        b$m_t_comment[i] <- paste0("No reported test, ", table)
      }
    }
    
    # Reported test exists
    else {
      # Special case: midp uses two values
      if (b$stat[i] == "midp") {
        if (!is.na(b$p_g[i])) {
          b$m_t[i] <- ifelse(b$p_match1[i] > b$p_g[i] & b$p_match2[i] > b$p_g[i], 
                             "Match", "No match")
        } else if (!is.na(b$p_l[i])) {
          b$m_t[i] <- ifelse(b$p_match1[i] < b$p_l[i] & b$p_match2[i] < b$p_l[i], 
                             "Match", "No match")
        }
        b$m_t_comment[i] <- "Midp exact&classic"
        
      } else {
        # Standard test
        if (!is.na(b$p_g[i])) {
          b$m_t[i] <- ifelse(b$p_match1[i] > b$p_g[i], "Match", "No match")
        } else if (!is.na(b$p_l[i])) {
          b$m_t[i] <- ifelse(b$p_match1[i] < b$p_l[i], "Match", "No match")
        }
        b$m_t_comment[i] <- b$stat[i]
      }
    }
  }
  
  # Capitalize first letter of comments
  b$m_t_comment <- ifelse(!is.na(b$m_t_comment), paste0(toupper(substr(b$m_t_comment, 1, 1)), 
                                                        substr(b$m_t_comment, 2, nchar(b$m_t_comment))), NA)
  
  #are there matching reported pvalues
  nm_ph <- c("p_chi", "p_fisher", "p_midp_exact", "p_midp_classic", "p_chi_yates", "p_lr", "p_cmh")
  names (nm_p) <- c("chisq", "fisher", "chisqc", "lr", "mh", "midp_exact","midp_classic")
  b$m_rep <- NA
  b$m_any <- NA
  
  for (i in seq_len(nrow(b))) {
    if (is.na(b$p_n[i])) next  # Skip if no reported numeric p-value
    
    #check all tests for exact match
    p_calc <- as.numeric(b[i, nm_ph])
    matches <- abs(b$p_n[i] - p_calc) < 1e-9
    
    if (any(matches, na.rm = TRUE)) {# Found at least one exact match
      match_idx <- which(matches)[1]
      b$m_any [i] <- names(nm_p)[nm_p == nm_ph[match_idx]]
    } else {(b$m_any [i] <- "No match")}
    
    if (!is.na(b$stat[i])) { # check specific test
      if (b$stat[i] == "midp") {
        b$m_rep [i] <- ifelse (abs(b$p_n[i] - b$p_midp_exact[i]) < 1e-9, "Midp_exact", 
                               ifelse (abs(b$p_n[i] - b$p_midp_classic[i]) < 1e-9, "Midp_classic", NA))
        
      } else {
        b$m_rep [i] <- ifelse (abs(b$p_n[i] - b[i, nm_p[b$stat[i]]]) < 1e-9, b$stat[i], NA)}
      b$m_rep [i] <- ifelse (is.na(b$m_rep [i]), "No match", b$m_rep [i])  
    }}
  
  b$m_rep <- ifelse(!is.na(b$m_rep), paste0(toupper(substr(b$m_rep, 1, 1)), 
                                            substr(b$m_rep, 2, nchar(b$m_rep))), NA)
  b$m_any <- ifelse(!is.na(b$m_any), paste0(toupper(substr(b$m_any, 1, 1)), 
                                            substr(b$m_any, 2, nchar(b$m_any))), NA)
  
  #calculate missing p-values
  b <- dplyr::left_join(b, a %>% dplyr::filter(calculated_level == "yes") %>% 
                          dplyr::select(study, var, calculated_level) %>% dplyr::distinct(),
                        by = c("study", "var"))
  
  b$change_made <- b$data_missing_match <- b$p_missing_match <-  NA
  
  # Function to calculate p-value based on test type
  calc_p <- function(mat, method) {
    p <- p1 <- p2 <- NA
    if (method == "chisq") {p <- chi (mat)$p
    } else if (method == "fisher") { p <- fisher_fast(mat)
    } else if (method == "chisqc") { p <-  chi(mat, correct = TRUE)$p
    } else if (method == "lr") { p <- lr(mat)
    } else if (method == "mh") { p <- cmh (mat)
    } else if (method == "midp") {
      p <- NA; pf <- fisher_fast(mat); p1 <- midp_exact(mat); p2 <- midp_classic(mat, pf)
    } else return("Incorrect statistic- see description for valid tests")
    return (list(p = p, p1 =p1, p2= p2))
  }
  
  # Function to get best matching p-value
  get_best_p <- function(result, target_p, digits) {
    p_values <- c(result$p, result$p1, result$p2)
    p_values <- p_values[!is.na(p_values)]
    
    if (length(p_values) == 0) return(list(p = NA, diff = Inf, matched = FALSE))
    
    p_values_rounded <- round(p_values, digits)
    diffs <- abs(p_values_rounded - target_p)
    best_idx <- which.min(diffs)
    
    return(list(p = p_values_rounded[best_idx], 
                diff = diffs[best_idx],
                matched = (diffs[best_idx] == 0)))
  }
  
  # Helper function to test a matrix adjustment
  test_adjustment <- function(test_x, change_description, current_best_diff, 
                              stat_method, reported_p, p_digits) {
    test_result <- calc_p(test_x, stat_method)
    test_best <- get_best_p(test_result, reported_p, p_digits)
    
    list(matrix = test_x, p = test_best$p, diff = test_best$diff, matched = test_best$matched,
         change = change_description, improved = test_best$diff < current_best_diff)}
  
  for (i in 1:nrow(b)) {
    if (b$m_rep[i] %in% "No match" & !is.na(b$stat[i]) & !is.na(b$p_n[i]) & b$calculated_level[i] %in% "yes" &
        ((b$levels_no [i]* b$group [i] >4 & startsWith(tolower(b$stat [i]), "chisq"))  | b$levels_no [i]* b$group [i] ==4)) {
      counts <- na.omit(as.numeric(b[i, n_cols]))
      x <- matrix(counts, byrow = TRUE, nrow = b$levels_no[i])
      
      reported_p <- b$p_n[i]
      stat_method <- tolower(b$stat[i])
      p_digits <- b$p_d[i]
      
      # Initial p-value and difference
      current_result <- calc_p(x, stat_method)
      current_best <- get_best_p(current_result, reported_p, p_digits)
      
      best_diff <- current_best$diff
      best_p <- current_best$p
      best_matrix <- x
      best_change <- "none"
      found_match <- FALSE
      
      n_groups <- ncol(x)
      
      for (delta in 1:min (x [1,])) {
        delta_best_diff <- Inf
        
        # Generate all combinations that sum to delta
        combs <- expand.grid(rep(list(0:delta), n_groups))
        combs <- combs[rowSums(combs) > 0, ]
        
        # Try each combination
        for (j in 1:nrow(combs)) {
          test_x <- x
          test_x[1, ] <- test_x[1, ] - as.numeric(combs[j, ])
          
          # Build description
          changes <- sapply(1:n_groups, function(g) {
            if (combs[j, g] > 0) paste0(combs[j, g], " from group ", g)
            else NULL })
          changes <- changes[!sapply(changes, is.null)]
          change_desc <- paste0("Removed ", paste(changes, collapse = ", "))
          
          result <- test_adjustment(test_x, change_desc, best_diff, 
                                    stat_method, reported_p, p_digits)
          
          delta_best_diff <- min(delta_best_diff, result$diff)
          
          if (result$improved) {
            best_diff <- result$diff
            best_p <- result$p
            best_matrix <- result$matrix
            best_change <- result$change
            
            if (result$matched) {
              found_match <- TRUE
              break}
          }
        }
        
        if (found_match) break
        
        # Stop if this delta level made things worse than baseline
        if (delta > 1 && delta_best_diff >= current_best$diff) {break}
      }
      
      if (best_change == "none") {b$change_made[i] <- "No match with missing data"}
      
      if (best_change != "none") {
        b$data_missing_match[i] <- paste(as.vector(t(best_matrix)), collapse = ", ")
        b$change_made[i] <- best_change
        b$p_missing_match[i] <- best_p}
    }
  }  
  
  b$p_d <- dplyr::coalesce(b$p_d, 3)
  formats <- sprintf("%%.%df", b$p_d)
  
  b <- b %>% dplyr::mutate(across(c(p_chi:p_midp_classic, p_missing_match),
                                  ~ sprintf(formats, .x)))
  
  #final dataset
  b <- b %>% dplyr::select(study, var, p, stat, m_any, m_rep, m_t, m_t_comment, p_chi:chi_warning,
                           calculated_level, p_missing_match, data_missing_match, change_made) %>%
    dplyr::left_join(a %>% dplyr::select(study, var, N_flag) %>% dplyr::distinct(), by = c("study", "var")) %>% 
    dplyr::rename(reported_pvalue= p, statistical_test= stat, any_match = m_any, match_reported_p = m_rep, 
                  threshold_match= m_t, threshold_test = m_t_comment) %>% 
    dplyr::arrange(match(study, unique(a$study)), match(var, unique(a$var)))
  
  return (list(b, max_digits))
} 


# shiny app ---------------------------------------------------------------


ui <- fluidPage(
  
  titlePanel("Check p-values for baseline categorical variables from a randomised controlled trial"),
  
  sidebarLayout(
    sidebarPanel(
      width= 3,
      h4("Data Structure"),
      numericInput("groups", "Number of Groups:", value = 2, min = 1, max =20),
      numericInput("rows", "Number of Rows:", value = 10, min = 1, max=1000),
      actionButton("configure", "Enter data", width = "100%"),
      hr(),
      
      conditionalPanel(condition = "input.configure > 0",
                       h4("Select Columns"),
                       tags$div(
                         tags$label("Required columns (always included):"),
                         tags$div(style = "margin-left: 20px; margin-bottom: 10px; color: #555;",
                                  tags$div("✓ Study"), tags$div("✓ Variable"),
                                  tags$div("✓ n"), tags$div("✓ N"), tags$div("✓ P-value"))),
                       
                       checkboxGroupInput("optional_columns", "Optional columns:",
                                          choices = list("Level" = "level", "Stat" = "stat", "Calculated_level" = "calculated_level"),
                                          selected = NULL), # change )), to ), 
                       
                       conditionalPanel(        # Column arrangement for n/N only
                         condition = "input.configure > 0",
                         hr(),
                         h4("n/N Column Pattern"),
                         radioButtons("nN_pattern",
                                      "Choose n and N arrangement:",
                                      choices = list("n1, N1, n2, N2, ..." = "alternating",
                                                     "n1, n2, ..., N1, N2, ..." = "grouped"),
                                      selected = "alternating")),
                       hr(),
                       
                       h4("Actions"),
                       actionButton("create_table", "Create Table", width = "100%"), br(), br(),
                       actionButton("analyze", "Run Analysis", width = "100%"), br(), br(),
                       downloadButton("download", "Download Results"),br(), br(),
                       downloadButton("download_data", "Download Data"),br(), br(),
                       actionButton("reset", "Reset App", width = "100%"),
                       
                       div(style = "background-color: #f8f9fa; padding: 10px; margin-top: 20px; border-radius: 5px; font-size: 12px;",
                           strong("Feedback"), br(),
                           "For bugs email ",
                           a("m.bolland@auckland.ac.nz", href = "mailto:m.bolland@auckland.ac.nz"),
                           br(),
                           "I'd also love to hear if it has been useful for you."))
    ),
    
    mainPanel(width = 9,
              # Main data entry table
              conditionalPanel(
                condition = "input.create_table > 0", hr(), 
                h3("Data Entry Table"),
                uiOutput("table_summary"), br(),
                rhandsontable::rHandsontableOutput("data_entry_table"),
                p(style = "color: #666; margin-top: 10px;",
                  "• Fill in all cells with your data. Delete or paste over the examples.", br(),
                  "• You can copy and paste from Excel or type directly into the cells.", br(), 
                  "• You can copy and paste multiple rows at once.", br(), 
                  "• Empty cells are allowed for levels, stat, and P only, otherwise an error will occur", br(),
                  "• Choose the order of the columns using the button.", br(),
                  "• n1 is the number in group1 for the level (eg 10 males in example 1), N1 is the total number in group 1 (eg 20) or the sum of n1 for all levels if there are missing data", br(),
                  "• Variable should be the same for each variable, level is optional and can describe the level eg variable gender, level male in example 1. If there are 2 levels (Male, Female), you can enter only 1 level as in example 2 and the second level is calculated automatically", br (),
                  "• Permitted stats are 'chisq' for standard chisquare, 'chisqc' for chisq + continuity correction, 'fisher' for fisher's exact, 'lr' for likelihood ratio, 'mh' for CMH test, 'midp' for midp tests", br(),
                  "• Permitted p-values are 'ns' for not-significant, '>x' for > any value, '<x' for < any value, or the reported number", br(),
                  "• The column 'calculated_levels' is for any calculated levels manually eg if male is x, then female is N minus x"),
                br(),
                # Display results dataframe
                conditionalPanel(
                  condition = "output.results_available",
                  h3("Analysis Results"),
                  tableOutput("results_table"),
                  br(),  # Add some spacing
                  # Add explanatory text
                  p(style = "color: #666; margin-top: 10px;",
                    "Notes:", br(),
                    "• reported_pvalue: p-value from above table", br(),
                    "• any_match: the matching statistical test if it matches any calculated p-value", br(),
                    "• match_reported_p: if reported p-value matches the calculated p-value using the reported statistical test", br(),
                    "• threshold match: if the p-value is a threshold, it indicates whether it matches the calculated p-value for the reported test or any p-value when no reported statistical test", br(),
                    "• threshold test: statistical test used to compare p-value with threshold", br(),
                    "• p_chi to p_midp: calculated p-values using various tests, midp_exact (eg epitools in R), midp_classic (eg midp in SAS) see code for calculations", br(),
                    "• chi_warning: if expected cell <5 for 2*2 or the relevant rule for >2*2 tables", br(),
                    "• p_missing_match: matching p_value from progressively removing a maximum of 1/group for each group then 2/group etc; data_missing_match: first matching set of n/N giving matching p_value; change_made: changes made to N to get matching set", br(),
                    "• N_flag: if sum of n does not equal sum of N for any variable, sum of n replaces N"))
              ))))


server <- function(input, output, session) {
  
  main_data <- reactiveVal()  # Reactive value to store main data
  values <- reactiveValues(analysis_results = NULL, max_digits = 3)
  
  # Generate column order based on selections
  get_column_order <- reactive({
    required_cols <- c("study", "variable", "n", "N", "p")
    optional_cols <- input$optional_columns
    selected_cols <- c(required_cols, optional_cols)
    n_groups <- input$groups
    
    col_list <- c()
    col_list <- c(col_list, "Study", "Variable")
    if ("level" %in% selected_cols) col_list <- c(col_list, "level")
    has_N <- "N" %in% selected_cols    # Add n and N columns based on pattern
    if (has_N) {
      if (input$nN_pattern == "alternating") {      # Both n and N selected - use pattern
        for (i in 1:n_groups) {
          col_list <- c(col_list, paste0("n", i), paste0("N", i))}  # n1, N1, n2, N2, ...
      } else {  # n1, n2, ..., N1, N2, ...
        for (i in 1:n_groups) { col_list <- c(col_list, paste0("n", i))}
        for (i in 1:n_groups) { col_list <- c(col_list, paste0("N", i))}
      }
    } else { # Only n columns (always required)
      for (i in 1:n_groups) {col_list <- c(col_list, paste0("n", i))}
    }
    
    col_list <- c(col_list, "p")
    if ("stat" %in% selected_cols) col_list <- c(col_list, "stat")
    if ("calculated_level" %in% selected_cols) col_list <- c(col_list, "calculated_level")
    return(col_list)
  })
  
  # Create main data entry table
  observe({
    trigger <- input$create_table
    req(trigger > 0)
    
    main_data(NULL)
    req(get_column_order())
    n_groups <- input$groups 
    n_rows <- input$rows
    
    columns <- get_column_order()
    df <- data.frame(matrix("", nrow = n_rows, ncol = length(columns)))
    colnames(df) <- columns
    df[1, c("Study", "Variable", "n1", "n2", "p", "N1", "N2")] <- c("Example1", "Gender", 10, 15, "0.19", 20, 20) 
    if (n_rows > 1) {df[2, c("Study", "Variable", "n1", "n2", "p", "N1", "N2")] <- c("Example1", "Gender", 10, 5, "0.19", 20 , 20)} 
    if (n_rows > 2) {df[3, c("Study", "Variable", "n1", "n2", "N1", "N2")] <- c("Example2", "Gender", 10, 5, 20, 20)}
    if("level" %in% columns) {df$level[1:n_rows] <- c("M", "F") [n_rows]}
    if("stat" %in% columns) {df$stat[1:n_rows] <- "fisher"}
    n_cols <- paste0("N", 1:n_groups)
    main_data(df)
  })
  
  # Render main data entry table
  output$data_entry_table <- rhandsontable::renderRHandsontable({
    req(main_data())
    rhandsontable::rhandsontable(main_data(), rowHeaders = 1:nrow(main_data()),
                                 height = 400, stretchH = "all")
  })
  
  # Update main data when table is edited
  observeEvent(input$data_entry_table, {
    main_data(rhandsontable::hot_to_r(input$data_entry_table))
  })
  
  observeEvent(input$analyze, {  # Clear previous results first
    
    values$analysis_results <- NULL
    values$analysis_data <- NULL
    
    if (!is.null(main_data())) { # Check for missing values
      table_data <- main_data()
      table_data[table_data == ""] <- NA
      
      tryCatch({
        a <- base_p_cat_format(as.data.frame(table_data))
        
        if (is.character(a) && length(a) == 1) { stop(a) } # Throw error with the message from the function
        
        b <- base_p_cat_calcs(a)
        values$analysis_results <- b[[1]]
        values$analysis_data <- a
        values$max_digits <- b[[2]]}, 
        
        error = function(e) {values$analysis_results <- data.frame(Error = paste("Error in analysis:", e$message))
        values$max_digits <- 3})
    } else {
      showNotification("Please create a table first", type = "warning")
    }
  })
  
  # Display the results dataframe
  output$results_table <- renderTable({values$analysis_results}, digits = function() values$max_digits)
  
  output$results_available <- reactive({!is.null(values$analysis_results)})
  outputOptions(output, "results_available", suspendWhenHidden = FALSE)
  
  # Download handler
  output$download <- downloadHandler(
    filename = function() {paste("analysis_results_", Sys.Date(), ".csv", sep = "")},
    content = function(file) {write.csv(values$analysis_results, file, row.names = FALSE, na = "" )})
  
  output$download_data <- downloadHandler(
    filename = function() {paste("analysis_data_", Sys.Date(), ".csv", sep = "")},
    content = function(file) {write.csv(values$analysis_data, file, row.names = FALSE, na = "")})
  
  observeEvent(input$reset, {session$reload()})
}

# Run the application 
shinyApp(ui = ui, server = server)



