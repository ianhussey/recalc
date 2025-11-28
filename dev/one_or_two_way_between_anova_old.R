#' ANOVA sanity-check from summary statistics (1- or 2-way)
#'
#' Computes possible ranges of F statistics from one- or two-way between-subjects
#' ANOVA using reported cell means, SDs and Ns, accounting for rounding in
#' the reported means and SDs. Returns a data frame of F values, p-values,
#' partial eta-squared and confidence intervals.
#'
#' This function is adapted from Nick Brown's \code{f_range()} and
#' the \pkg{rpsychi} package, with modifications to work without external
#' dependencies and to return a tidy data frame.
#'
#' @param m Matrix of means (for 2-way ANOVA) or numeric vector (for 1-way).
#' @param sd Matrix of standard deviations (same shape as \code{m}) or numeric vector.
#' @param n Matrix or vector of sample sizes. If a vector, it is recycled to the
#'   shape of \code{m}.
#' @param show.t Logical. If \code{TRUE}, statistics are converted from F to t
#'   by taking square roots.
#' @param dp.p Integer. Number of decimal places assumed in the reported
#'   means and SDs. If \code{-1} (default), this is inferred from \code{m} and
#'   \code{sd}.
#' @param labels Optional character vector of labels for the omnibus tests.
#'   For 2-way ANOVA, the defaults are \code{c("col F", "row F", "inter F")};
#'   for 1-way ANOVA, \code{"F"} or \code{"t"} depending on \code{show.t}.
#' @param peta2_sig_level Significance level used when computing confidence
#'   intervals for partial eta-squared.
#' @param F_rounding_set Character vector of rounding methods applied to the
#'   reported/min/max F statistics. Allowed values are
#'   \code{"half_up"}, \code{"half_down"}, \code{"bankers"}, and \code{"trunc"}.
#'
#' @return A data frame with one row per combination of effect (row/col/interaction
#'   or 1-way), F-type (reported, min, max), and rounding method, containing:
#'   \code{F}, \code{df1}, \code{df2}, \code{p}, \code{peta2}, and
#'   confidence limits for \code{peta2}.
#'
#' @importFrom roundwork round_up round_down round_trunc
#'
#' @examples
#' m  <- matrix(c(2.93, 1.98,
#'                2.71, 1.84),
#'              ncol = 2)
#' sd <- matrix(c(1.10, 0.95,
#'                1.05, 0.98),
#'              ncol = 2)
#' n  <- matrix(c(40, 20,
#'                35, 10),
#'              ncol = 2)
#'
#' res <- one_or_two_way_between_anova(m = m, sd = sd, n = n)
#' head(res)
one_or_two_way_between_anova <- function(m,
                                         sd,
                                         n,
                                         show.t = FALSE,
                                         dp.p = -1,
                                         labels = c(),
                                         peta2_sig_level = 0.10,
                                         F_rounding_set = c("half_up",
                                                            "half_down",
                                                            "bankers",
                                                            "trunc")) {
  # validate rounding methods (parallel to independent_t_test)
  allowed_F_rounding <- c("half_up", "half_down", "bankers", "trunc")
  bad_round <- setdiff(F_rounding_set, allowed_F_rounding)
  if (length(bad_round) > 0L) {
    stop(
      "F_rounding_set contains invalid values: ",
      paste0(bad_round, collapse = ", "),
      "\nAllowed: ", paste0(allowed_F_rounding, collapse = ", ")
    )
  }
  F_rounding_set <- unique(F_rounding_set)
  
  m.ok <- m
  
  if (inherits(m.ok, "matrix")) {
    func <- ind_twoway_second
    useF <- c(3, 2, 4)
    default_labels <- c("col F", "row F", "inter F")
  } else {
    m.ok <- matrix(m)
    func <- ind_twoway_second
    useF <- 1
    default_labels <- c("F")
    if (show.t) {
      default_labels <- c("t")
    }
  }
  
  # determine decimal places if not specified
  dp <- dp.p
  if (dp.p == -1) {
    dp <- 0
    numbers <- c(m, sd)
    for (i in numbers) {
      if (i != roundwork::round_up(i, 0)) {
        dp <- max(dp, 1)
        j <- i * 10
        if (j != roundwork::round_up(j, 0)) {
          dp <- max(dp, 2)
        }
      }
    }
  }
  
  if (length(labels) == 0) {
    labels <- default_labels
  }
  
  # nominal ANOVA table and F values (no rounding error)
  anova_tab <- func(m = m.ok, sd = sd, n = n)$anova.table
  f.nom     <- anova_tab$F
  
  # rounding error bounds on inputs
  delta   <- (0.1^dp) / 2
  tiny.sd <- delta / 100
  
  sd.hi <- pmax(sd - delta, tiny.sd)
  sd.lo <- pmax(sd + delta, tiny.sd)
  
  f.hi <- f.nom
  f.lo <- f.nom
  
  # generate all combinations of +/- delta on the means
  l <- length(m.ok)
  rawcomb <- utils::combn(rep(c(-delta, delta), l), l)
  comb    <- rawcomb[, !duplicated(t(rawcomb))]
  
  for (i in seq_len(ncol(comb))) {
    m.adj <- m.ok + comb[, i]
    
    if (all((abs(m.adj - m.adj[1])) < 1e-14)) {
      adj.f.hi <- 0
      adj.f.lo <- 0
    } else {
      adj_tab_h <- func(m = m.adj, sd = sd.hi, n = n)$anova.table
      adj_tab_l <- func(m = m.adj, sd = sd.lo, n = n)$anova.table
      adj.f.hi  <- adj_tab_h$F
      adj.f.lo  <- adj_tab_l$F
    }
    
    f.hi <- pmax(f.hi, adj.f.hi)
    f.lo <- pmin(f.lo, adj.f.lo)
  }
  
  if (show.t) {
    f.nom <- sqrt(f.nom)
    f.hi  <- sqrt(f.hi)
    f.lo  <- sqrt(f.lo)
  }
  
  # result template
  result_df <- data.frame(
    label      = character(),
    F_type     = character(),
    F_rounding = character(),
    F          = numeric(),
    df1        = numeric(),
    df2        = numeric(),
    p          = numeric(),
    peta2      = numeric(),
    stringsAsFactors = FALSE
  )
  
  # df_effect and df_error for each omnibus test
  df_effect_vec <- anova_tab$df[useF]
  df_error      <- anova_tab$df[5L]  # "Within"
  
  # decimal places used when rounding F (keep previous behaviour)
  fdp <- 3
  
  # loop over effects (row, col, interaction or 1-way)
  for (i in seq_along(useF)) {
    j         <- useF[i]
    df_effect <- df_effect_vec[i]
    
    nominal_p <- 1 - stats::pf(f.nom[j], df_effect, df_error)
    min_p     <- 1 - stats::pf(f.lo[j],  df_effect, df_error)
    max_p     <- 1 - stats::pf(f.hi[j],  df_effect, df_error)
    
    peta2_nom <- (f.nom[j] * df_effect) /
      (f.nom[j] * df_effect + df_error)
    peta2_min <- (f.lo[j] * df_effect) /
      (f.lo[j] * df_effect + df_error)
    peta2_max <- (f.hi[j] * df_effect) /
      (f.hi[j] * df_effect + df_error)
    
    # iterate over F rounding schemes (parallel to d_rounding_set in independent_t_test)
    for (F_rounding in F_rounding_set) {
      F_round_fun <- switch(
        F_rounding,
        "half_up"   = function(x) roundwork::round_up(x,   fdp),
        "half_down" = function(x) roundwork::round_down(x, fdp),
        "bankers"   = function(x) round(x, fdp),
        "trunc"     = function(x) roundwork::round_trunc(x, fdp)
      )
      
      result_df <- rbind(
        result_df,
        data.frame(
          label      = labels[i],
          F_type     = "reported",
          F_rounding = F_rounding,
          F          = F_round_fun(f.nom[j]),
          df1        = df_effect,
          df2        = df_error,
          p          = nominal_p,
          peta2      = peta2_nom
        ),
        data.frame(
          label      = labels[i],
          F_type     = "min",
          F_rounding = F_rounding,
          F          = F_round_fun(f.lo[j]),
          df1        = df_effect,
          df2        = df_error,
          p          = min_p,
          peta2      = peta2_min
        ),
        data.frame(
          label      = labels[i],
          F_type     = "max",
          F_rounding = F_rounding,
          F          = F_round_fun(f.hi[j]),
          df1        = df_effect,
          df2        = df_error,
          p          = max_p,
          peta2      = peta2_max
        )
      )
    }
  }
  
  # add CIs for partial eta-squared
  ci_df <- result_df |>
    dplyr::rowwise() |>
    dplyr::do(petasq_confint(
      petasq    = .$peta2,
      df1       = .$df1,
      df2       = .$df2,
      sig_level = peta2_sig_level
    )) |>
    dplyr::ungroup() |>
    dplyr::select(
      peta2_ci_lower = partial_etasq_lower,
      peta2_ci_upper = partial_etasq_upper
    )
  
  dplyr::bind_cols(result_df, ci_df)
}

#' Sample variance (population denominator)
#'
#' Computes the variance of a numeric vector using denominator \eqn{n}
#' (rather than \eqn{n - 1}), i.e. the mean squared deviation from the mean.
#'
#' @param x Numeric vector.
#' @param na_rm Logical. If \code{TRUE}, remove \code{NA}s before computing.
#'
#' @return A single numeric value giving the variance.
#' @keywords internal
svar <- function(x, na_rm = TRUE) {
  if (na_rm) {
    x <- stats::na.omit(x)
  }
  1 / length(x) * sum((x - mean(x))^2)
}

#' Convert sum of squared deviations to standard deviation
#'
#' Converts a sum of squared deviations (SSD) to a standard deviation
#' assuming an unbiased estimator.
#'
#' @param n Sample size.
#' @param ssd Sum of squared deviations.
#'
#' @return Standard deviation.
#' @keywords internal
ssd2sd <- function(n, ssd) {
  sqrt((n / (n - 1)) * ssd^2)
}

#' Noncentrality parameter for the F distribution
#'
#' Finds the noncentrality parameter of an F distribution that yields
#' a given cumulative probability at a specified F value.
#'
#' @param x Observed F statistic.
#' @param df1 Numerator degrees of freedom.
#' @param df2 Denominator degrees of freedom.
#' @param prob Target cumulative probability.
#' @param interval Numeric vector of length 2 giving the search interval.
#' @param my_tol Tolerance passed to \code{\link[stats]{uniroot}}.
#'
#' @return Estimated noncentrality parameter (numeric scalar).
#' @keywords internal
F_test_noncentrality <- function(x, df1, df2,
                                 prob,
                                 interval = c(0, 10000),
                                 my_tol = 1e-06) {
  f_root <- function(ncp) stats::pf(x, df1, df2, ncp) - prob
  stats::uniroot(f_root, interval, tol = my_tol)$root
}

#' Approximate power for an F-test given Cohen's f
#'
#' Computes the power of an F-test in terms of Cohen's f effect size
#' using the approximation from \emph{rpsychi}.
#'
#' @param df1 Numerator degrees of freedom.
#' @param df2 Denominator degrees of freedom.
#' @param delta Cohen's f (or vector of f values).
#' @param sig.level Significance level (default 0.05).
#'
#' @return Numeric vector of power values (same length as \code{delta}).
#' @keywords internal
power_f2 <- function(df1, df2, delta, sig.level = 0.05) {
  n <- df2 / (df1 + 1) + 1
  fc <- stats::qf(p = sig.level,
                  df1 = df1,
                  df2 = (df1 + 1) * (n - 1),
                  lower.tail = FALSE)
  lamda <- (delta^2) * (n * (df1 + 1))
  v <- (df1 + 1) * (n - 1)
  z1b <- (sqrt(2 * (df1 + lamda) - ((df1 + 2 * lamda) / (df1 + lamda))) -
            sqrt((2 * v - 1) * ((df1 * fc) / v))) /
    sqrt(((df1 * fc) / v) + ((df1 + 2 * lamda) / (df1 + lamda)))
  stats::pnorm(z1b)
}

#' Partial eta-squared and confidence intervals
#'
#' Computes partial eta-squared and its confidence interval using
#' the noncentral F distribution.
#'
#' @param petasq Vector of partial eta-squared values.
#' @param df1 Vector of numerator degrees of freedom (same length as \code{petasq}).
#' @param df2 Denominator degrees of freedom (scalar).
#' @param sig_level Significance level for the confidence interval
#'   (default 0.10, giving 90\% CIs).
#'
#' @return A data frame with columns
#'   \code{partial_etasq}, \code{partial_etasq_lower}, \code{partial_etasq_upper}.
#' @keywords internal
petasq_confint <- function(petasq, df1, df2, sig_level = 0.10) {
  f2 <- petasq / (1 - petasq)
  f.value <- f2 * (df2 / df1)
  iter <- length(petasq)
  delta_lower <- delta_upper <- numeric(iter)
  
  for (i in seq_len(iter)) {
    delta_lower[i] <- try(
      F_test_noncentrality(
        f.value[i], df1[i], df2,
        prob = 1 - sig_level / 2
      ),
      silent = TRUE
    )
    delta_upper[i] <- try(
      F_test_noncentrality(
        f.value[i], df1[i], df2,
        prob = sig_level / 2
      ),
      silent = TRUE
    )
  }
  
  cond1 <- is.character(delta_lower)
  cond2 <- is.character(delta_upper)
  
  if (cond1) {
    delta_lower[grep("Error", delta_lower)] <- 0
    delta_lower <- as.numeric(delta_lower)
  }
  
  if (cond2) {
    delta_upper[grep("Error", delta_upper)] <- 0
    delta_upper <- as.numeric(delta_upper)
  }
  
  lower_petasq <- delta_lower / (delta_lower + df1 + df2 + 1)
  upper_petasq <- delta_upper / (delta_upper + df1 + df2 + 1)
  
  data.frame(
    partial_etasq        = petasq,
    partial_etasq_lower  = lower_petasq,
    partial_etasq_upper  = upper_petasq
  )
}

#' Two-way between-subjects ANOVA from summary statistics
#'
#' Computes a two-way between-subjects ANOVA from cell means, standard
#' deviations, and sample sizes. The implementation is adapted from the
#' \pkg{rpsychi} package and Nick Brown's \code{f_range()} code.
#'
#' @param m Matrix of cell means (rows = levels of factor A, columns = levels of factor B).
#' @param sd Matrix of cell standard deviations (same shape as \code{m}).
#' @param n Matrix (or vector) of cell sample sizes; recycled to the
#'   shape of \code{m} if necessary.
#' @param unbiased Logical. If \code{FALSE}, treat \code{sd} as sum of
#'   squared deviations and convert to unbiased SD. Default \code{TRUE}.
#' @param sig_level Significance level for partial eta-squared confidence intervals
#'   (default 0.05).
#' @param digits Number of decimal places for rounding in the returned
#'   objects (currently not applied; retained for compatibility).
#'
#' @return A list with components:
#' \describe{
#'   \item{anova.table}{\code{data.frame} with SS, df, MS, F for
#'         total between, row, column, interaction, within, and total.}
#'   \item{omnibus.es}{\code{data.frame} with partial eta-squared and its
#'         confidence intervals for row, column, and interaction effects.}
#'   \item{power}{\code{matrix} with approximate power for small, medium,
#'         and large effect sizes for each effect.}
#' }
#' @export
ind_twoway_second <- function(m,
                              sd,
                              n,
                              unbiased = TRUE,
                              sig_level = 0.05,
                              digits = 3) {
  
  # ensure n is a matrix
  if (is.vector(n)) {
    n <- matrix(n, ncol = ncol(m), nrow = nrow(m))
  }
  
  # harmonic N
  if (nlevels(as.factor(n)) == 1) {
    Nh <- sum(n)
  } else {
    Nh <- (prod(dim(m))^2) / sum(1 / n)
  }
  
  # convert SSD to SD if requested
  if (!unbiased) {
    sd <- ssd2sd(n, sd)
  }
  
  dfw <- (sum(n) - prod(dim(m)))
  MSw <- sum((n - 1) * sd^2) / dfw
  SSb <- Nh * svar(as.vector(m))
  SSrow <- Nh * svar(rowSums(m) / ncol(m))
  SScol <- Nh * svar(colSums(m) / nrow(m))
  SSint <- SSb - SSrow - SScol
  SSw <- MSw * dfw
  SSt <- SSb + SSw
  
  dfrow <- nrow(m) - 1
  dfcol <- ncol(m) - 1
  dfint <- dfrow * dfcol
  dfb <- dfrow + dfcol + dfint
  
  MSrow <- SSrow / dfrow
  MScol <- SScol / dfcol
  MSint <- SSint / dfint
  
  frow <- MSrow / MSw
  fcol <- MScol / MSw
  fint <- MSint / MSw
  
  MSb <- SSb / dfb
  fb <- MSb / MSw
  
  # p-values
  p_row <- stats::pf(q = frow, df1 = dfrow, df2 = dfw, lower.tail = FALSE)
  p_col <- stats::pf(q = fcol, df1 = dfcol, df2 = dfw, lower.tail = FALSE)
  p_int <- stats::pf(q = fint, df1 = dfint, df2 = dfw, lower.tail = FALSE)
  p_b   <- stats::pf(q = fb,   df1 = dfb,  df2 = dfw, lower.tail = FALSE)
  
  # ANOVA table
  anova.table <- data.frame(
    SS = c(SSb,   SSrow, SScol, SSint, SSw,  SSt),
    df = c(dfb,   dfrow, dfcol, dfint, dfw,  dfb + dfw),
    MS = c(MSb,   MSrow, MScol, MSint, MSw,  NA_real_),
    F  = c(fb,    frow,  fcol,  fint,  NA,   NA),
    row.names = c("Between",
                  "Between (row)",
                  "Between (col)",
                  "Between (row * col)",
                  "Within",
                  "Total")
  )
  class(anova.table) <- c("anova", "data.frame")
  
  # partial eta-squared
  petasq_row <- SSrow / (SSrow + SSw)
  petasq_col <- SScol / (SScol + SSw)
  petasq_int <- SSint / (SSint + SSw)
  
  omnibus.es <- petasq_confint(
    petasq   = c(petasq_row, petasq_col, petasq_int),
    df1      = c(dfrow, dfcol, dfint),
    df2      = dfw,
    sig_level = sig_level
  )
  rownames(omnibus.es) <- c("Between (row)",
                            "Between (col)",
                            "Between (row * col)")
  
  # power for small/medium/large f
  c_delta <- c(0.1, 0.25, 0.4)
  criterion_power <- rbind(
    power_f2(sig.level = sig_level, df1 = dfrow, df2 = dfw, delta = c_delta),
    power_f2(sig.level = sig_level, df1 = dfcol, df2 = dfw, delta = c_delta),
    power_f2(sig.level = sig_level, df1 = dfint, df2 = dfw, delta = c_delta)
  )
  colnames(criterion_power) <- c("small", "medium", "large")
  rownames(criterion_power) <- rownames(omnibus.es)
  
  list(
    anova.table = anova.table,
    omnibus.es  = omnibus.es,
    power       = criterion_power
  )
}

#' Round numeric to character, preserving trailing zeros
#'
#' Rounds a numeric vector to a given number of decimal places and returns
#' a character vector with trailing zeros retained.
#'
#' @param x Numeric vector.
#' @param digits Number of decimal places.
#'
#' @return Character vector.
#' @export
round_half_up_to_char <- function(x, digits = 2) {
  rounded <- round(x, digits)
  formatC(rounded, format = "f", digits = digits)
}

#' #' ANOVA sanity-check from summary statistics (1- or 2-way)
#' #'
#' #' Computes possible ranges of F statistics from one- or two-way between-subjects
#' #' ANOVA using reported cell means, SDs and Ns, accounting for rounding in
#' #' the reported means and SDs. Returns a data frame of F values, p-values,
#' #' partial eta-squared and confidence intervals.
#' #'
#' #' This function is adapted from Nick Brown's \code{f_range()} and
#' #' the \pkg{rpsychi} package, with modifications to work without external
#' #' dependencies and to return a tidy data frame.
#' #'
#' #' @param m Matrix of means (for 2-way ANOVA) or numeric vector (for 1-way).
#' #' @param sd Matrix of standard deviations (same shape as \code{m}) or numeric vector.
#' #' @param n Matrix or vector of sample sizes. If a vector, it is recycled to the
#' #'   shape of \code{m}.
#' #' @param show.t Logical. If \code{TRUE}, statistics are converted from F to t
#' #'   by taking square roots.
#' #' @param dp.p Integer. Number of decimal places assumed in the reported
#' #'   means and SDs. If \code{-1} (default), this is inferred from \code{m} and
#' #'   \code{sd}.
#' #' @param labels Optional character vector of labels for the effects. If empty,
#' #'   reasonable defaults are constructed.
#' #' @param peta2_sig_level Significance level for partial eta-squared confidence
#' #'   intervals (default 0.10, giving 90\% CIs).
#' #'   
#' #' @importFrom roundwork round_up round_down
#' #' @importFrom dplyr bind_cols rowwise do ungroup select
#' #'
#' #' @return A data frame with one row per F variant and columns:
#' #'   \itemize{
#' #'     \item \code{label} Effect label (e.g., "Between (row)").
#' #'     \item \code{F_type} One of "reported", "min", "max".
#' #'     \item \code{F} F statistic (or t if \code{show.t = TRUE}).
#' #'     \item \code{df1}, \code{df2} Degrees of freedom.
#' #'     \item \code{p} p-value.
#' #'     \item \code{peta2} Partial eta-squared.
#' #'     \item \code{peta2_ci_lower}, \code{peta2_ci_upper} CI bounds for partial eta-squared.
#' #'   }
#' #' @export
#' #'
#' #' @examples
#' #' m <- matrix(c(5.00, 2.69,
#' #'               4.83, 5.54),
#' #'             ncol = 2)
#' #' sd <- matrix(c(2.99, 2.57,
#' #'                2.71, 1.84),
#' #'              ncol = 2)
#' #' n  <- matrix(c(40, 20,
#' #'                35, 10),
#' #'              ncol = 2)
#' #'
#' #' res <- one_or_two_way_between_anova(m = m, sd = sd, n = n)
#' #' head(res)
#' one_or_two_way_between_anova <- function(m,
#'                                      sd,
#'                                      n,
#'                                      show.t = FALSE,
#'                                      dp.p = -1,
#'                                      labels = c(),
#'                                      peta2_sig_level = 0.10) {
#'   # note: depend explicitly on janitor and dplyr
#'   m.ok <- m
#'   
#'   if (inherits(m.ok, "matrix")) {
#'     func <- ind_twoway_second
#'     useF <- c(3, 2, 4)
#'     default_labels <- c("col F", "row F", "inter F")
#'   } else {
#'     m.ok <- matrix(m)
#'     func <- ind_twoway_second
#'     useF <- 1
#'     default_labels <- c("F")
#'     if (show.t) {
#'       default_labels <- c("t")
#'     }
#'   }
#'   
#'   # determine decimal places if not specified
#'   dp <- dp.p
#'   if (dp.p == -1) {
#'     dp <- 0
#'     numbers <- c(m, sd)
#'     for (i in numbers) {
#'       if (i != janitor::round_half_up(i, 0)) {
#'         dp <- max(dp, 1)
#'         j <- i * 10
#'         if (j != janitor::round_half_up(j, 0)) {
#'           dp <- max(dp, 2)
#'         }
#'       }
#'     }
#'   }
#'   
#'   if (length(labels) == 0) {
#'     labels <- default_labels
#'   }
#'   
#'   # nominal F values (no rounding error)
#'   f.nom <- func(m = m.ok, sd = sd, n = n)$anova.table$F
#'   
#'   # rounding error bounds
#'   delta <- (0.1^dp) / 2
#'   tiny.sd <- delta / 100
#'   
#'   sd.hi <- pmax(sd - delta, tiny.sd)
#'   sd.lo <- pmax(sd + delta, tiny.sd)
#'   
#'   f.hi <- f.nom
#'   f.lo <- f.nom
#'   
#'   # generate all combinations of +/- delta on the means
#'   l <- length(m.ok)
#'   rawcomb <- utils::combn(rep(c(-delta, delta), l), l)
#'   comb <- rawcomb[, !duplicated(t(rawcomb))]
#'   
#'   for (i in seq_len(ncol(comb))) {
#'     m.adj <- m.ok + comb[, i]
#'     
#'     if (all((abs(m.adj - m.adj[1])) < 1e-14)) {
#'       adj.f.hi <- 0
#'       adj.f.lo <- 0
#'     } else {
#'       adj.f.hi <- func(m = m.adj, sd = sd.hi, n = n)$anova.table$F
#'       adj.f.lo <- func(m = m.adj, sd = sd.lo, n = n)$anova.table$F
#'     }
#'     
#'     f.hi <- pmax(f.hi, adj.f.hi)
#'     f.lo <- pmin(f.lo, adj.f.lo)
#'   }
#'   
#'   if (show.t) {
#'     f.nom <- sqrt(f.nom)
#'     f.hi  <- sqrt(f.hi)
#'     f.lo  <- sqrt(f.lo)
#'   }
#'   
#'   # result template
#'   result_df <- data.frame(
#'     label = character(),
#'     F_type = character(),
#'     F = numeric(),
#'     df1 = numeric(),
#'     df2 = numeric(),
#'     p = numeric(),
#'     peta2 = numeric(),
#'     stringsAsFactors = FALSE
#'   )
#'   
#'   df_effect <- nrow(m.ok) - 1
#'   df_error  <- sum(n) - nrow(m.ok)
#'   
#'   fdp <- 3
#'   
#'   for (i in seq_along(useF)) {
#'     j <- useF[i]
#'     
#'     nominal_p <- 1 - stats::pf(f.nom[j], df_effect, df_error)
#'     min_p     <- 1 - stats::pf(f.lo[j],  df_effect, df_error)
#'     max_p     <- 1 - stats::pf(f.hi[j],  df_effect, df_error)
#'     
#'     peta2_nom <- (f.nom[j] * df_effect) /
#'       (f.nom[j] * df_effect + df_error)
#'     peta2_min <- (f.lo[j] * df_effect) /
#'       (f.lo[j] * df_effect + df_error)
#'     peta2_max <- (f.hi[j] * df_effect) /
#'       (f.hi[j] * df_effect + df_error)
#'     
#'     result_df <- rbind(
#'       result_df,
#'       data.frame(
#'         label = labels[i],
#'         F_type = "reported",
#'         F = janitor::round_half_up(f.nom[j], fdp),
#'         df1 = df_effect,
#'         df2 = df_error,
#'         p = nominal_p,
#'         peta2 = peta2_nom
#'       ),
#'       data.frame(
#'         label = labels[i],
#'         F_type = "min",
#'         F = janitor::round_half_up(f.lo[j], fdp),
#'         df1 = df_effect,
#'         df2 = df_error,
#'         p = min_p,
#'         peta2 = peta2_min
#'       ),
#'       data.frame(
#'         label = labels[i],
#'         F_type = "max",
#'         F = janitor::round_half_up(f.hi[j], fdp),
#'         df1 = df_effect,
#'         df2 = df_error,
#'         p = max_p,
#'         peta2 = peta2_max
#'       )
#'     )
#'   }
#'   
#'   # add partial eta-squared CIs
#'   ci_df <- result_df |>
#'     dplyr::rowwise() |>
#'     dplyr::do(petasq_confint(
#'       petasq   = .$peta2,
#'       df1      = .$df1,
#'       df2      = .$df2,
#'       sig_level = peta2_sig_level
#'     )) |>
#'     dplyr::ungroup() |>
#'     dplyr::select(
#'       peta2_ci_lower = partial_etasq_lower,
#'       peta2_ci_upper = partial_etasq_upper
#'     )
#'   
#'   dplyr::bind_cols(result_df, ci_df)
#' }