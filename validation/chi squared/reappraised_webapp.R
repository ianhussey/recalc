# reappraised_webapp.R
# -----------------------------------------------------------------------------
# Per-table categorical p-value tests extracted from the reappraised Shiny app
# (apps/app_dichotomous.R, function base_p_cat_calcs()), MIT-licensed,
# Bolland et al. The web app is the most capable reappraised tool: it computes
# SEVEN tests per table and flags a reported p as reproducible if it matches
# ANY of them at the reported precision (its `m_any` / "match any" rule). That
# is a gap-aware union of discrete test values -- directly analogous to recalc's
# union, but over the app's seven tests (it lacks only recalc's central Fisher
# and unpooled Wald z).
#
# The seven tests, mapped to their recalc equivalents:
#   p_chi          -> pearson            (Pearson, uncorrected)
#   p_chi_yates    -> yates              (Yates continuity correction)
#   p_fisher       -> fisher             (Fisher exact, two-sided minlike)
#   p_lr           -> likelihood_ratio   (G^2)
#   p_cmh          -> n1_chisq           (CMH general association = (N-1)/N * Pearson for 2x2)
#   p_midp_exact   -> fisher_midp        (central mid-p, epitools convention)
#   p_midp_classic -> fisher_midp_sas    (SAS PROC FREQ mid-p)
#
# The functions are lifted with only the data-frame plumbing and the
# expected-cell warning removed; the numerics are unchanged from the app. `x`
# is a 2x2 integer matrix (rows = groups, columns = level 1 / level 2).
# -----------------------------------------------------------------------------

# Chi-square (uncorrected or Yates-corrected) -- verbatim numerics from app `chi()`
.webapp_chi <- function(x, correct = FALSE) {
  if (nrow(x) == 2L && ncol(x) == 2L) {
    # fast 2x2
    a <- x[1, 1]
    b <- x[1, 2]
    c <- x[2, 1]
    d <- x[2, 2]
    n <- a + b + c + d
    if (correct) {
      chi_stat <- n *
        (max(abs(a * d - b * c) - n / 2, 0))^2 /
        ((a + b) * (c + d) * (a + c) * (b + d))
    } else {
      chi_stat <- n *
        (a * d - b * c)^2 /
        ((a + b) * (c + d) * (a + c) * (b + d))
    }
    return(pchisq(chi_stat, 1L, lower.tail = FALSE))
  } else {
    # general r x c
    nr <- nrow(x)
    nc <- ncol(x)
    row <- .rowSums(x, nr, nc)
    col <- .colSums(x, nr, nc)
    n <- sum(row)
    exp <- tcrossprod(row, col) / n
    if (correct && nr == 2L && nc == 2L) {
      chi_stat <- sum((pmax(abs(x - exp) - 0.5, 0))^2 / exp)
    } else {
      diff <- x - exp
      chi_stat <- sum(diff * diff / exp)
    }
    df <- (nr - 1L) * (nc - 1L)
    return(pchisq(chi_stat, df, lower.tail = FALSE))
  }
}

# Fisher's exact -- verbatim from app `fisher_fast()`. 2x2 path is exact
# (two-sided minlike); larger tables use stats' Monte-Carlo C routine.
.webapp_fisher_fast <- function(x, alternative = "two.sided", B = 10000) {
  nr <- nrow(x)
  nc <- ncol(x)
  if (nr == 2L && nc == 2L) {
    a <- x[1, 1]
    b <- x[1, 2]
    c <- x[2, 1]
    d <- x[2, 2]
    m <- a + c
    n <- b + d
    k <- a + b
    if (alternative == "two.sided") {
      p_obs <- dhyper(a, m, n, k)
      lo <- max(0L, k - n)
      hi <- min(k, m)
      support <- lo:hi
      p_all <- dhyper(support, m, n, k)
      p_value <- sum(p_all[p_all <= p_obs * (1 + 1e-7)])
    } else if (alternative == "less") {
      p_value <- phyper(a, m, n, k)
    } else if (alternative == "greater") {
      p_value <- phyper(a - 1, m, n, k, lower.tail = FALSE)
    }
    return(p_value)
  }
  sr <- rowSums(x)
  sc <- colSums(x)
  x <- x[sr > 0, sc > 0, drop = FALSE]
  STATISTIC <- -sum(lfactorial(x))
  tmp <- .Call(stats:::C_Fisher_sim, rowSums(x), colSums(x), as.integer(B))
  almost.1 <- 1 + 64 * .Machine$double.eps
  p_value <- (1 + sum(tmp <= STATISTIC / almost.1)) / (B + 1)
  max(0, min(1, p_value))
}

# Likelihood-ratio G^2 -- verbatim from app `lr()`
.webapp_lr <- function(obs) {
  row <- rowSums(obs)
  col <- colSums(obs)
  n <- sum(obs)
  exp <- outer(row, col) / n
  g_stat <- 2 * sum(obs[obs > 0] * log(obs[obs > 0] / exp[obs > 0]))
  df <- (nrow(obs) - 1) * (ncol(obs) - 1)
  pchisq(g_stat, df, lower.tail = FALSE)
}

# CMH general-association test -- verbatim from app `cmh()`
# (for a 2x2 table this equals (N-1)/N * Pearson = recalc's n1_chisq)
.webapp_cmh <- function(x) {
  R <- nrow(x)
  C <- ncol(x)
  nt <- sum(x)
  pr <- rowSums(x) / nt
  pc <- colSums(x) / nt
  m <- as.vector(nt * outer(pr, pc))
  n <- as.vector(x)
  V1 <- diag(pr) - pr %*% t(pr)
  V2 <- diag(pc) - pc %*% t(pc)
  V <- (nt^2 / (nt - 1)) * kronecker(V2, V1)
  A_row <- cbind(diag(R - 1), rep(0, R - 1))
  A_col <- cbind(diag(C - 1), rep(0, C - 1))
  A <- kronecker(A_col, A_row)
  AVA <- A %*% V %*% t(A)
  Q <- as.numeric(t(n - m) %*% t(A) %*% solve(AVA) %*% A %*% (n - m))
  df <- (R - 1) * (C - 1)
  pchisq(Q, df, lower.tail = FALSE)
}

# Central mid-p (epitools convention) -- verbatim from app `midp_exact()`
.webapp_midp_exact <- function(x) {
  lteq <- .webapp_fisher_fast(x, alternative = "less")
  gteq <- .webapp_fisher_fast(x, alternative = "greater")
  pval1 <- 0.5 * (lteq - gteq + 1)
  one_sided <- min(pval1, 1 - pval1)
  2 * one_sided
}

# SAS PROC FREQ mid-p -- verbatim from app `midp_classic()`
.webapp_midp_classic <- function(x, p_fisher) {
  a <- x[1, 1]
  b <- x[1, 2]
  c <- x[2, 1]
  d <- x[2, 2]
  n <- sum(x)
  log_prob <- (lgamma(a + b + 1) +
    lgamma(c + d + 1) +
    lgamma(a + c + 1) +
    lgamma(b + d + 1) -
    lgamma(n + 1) -
    lgamma(a + 1) -
    lgamma(b + 1) -
    lgamma(c + 1) -
    lgamma(d + 1))
  prob <- exp(log_prob)
  p_fisher - 0.5 * prob
}

# -----------------------------------------------------------------------------
# Wrapper: the seven web-app test p-values for a 2x2 table, named by their
# recalc-equivalent method. The web app's reproducibility decision is a
# match-any (union) over these seven values at the reported precision.
# -----------------------------------------------------------------------------
reappraised_webapp_pvalues <- function(tab) {
  if (!(nrow(tab) == 2L && ncol(tab) == 2L)) {
    stop(
      "reappraised_webapp_pvalues() is defined for 2x2 tables (the app's mid-p path)."
    )
  }
  pf <- .webapp_fisher_fast(tab)
  c(
    pearson = as.numeric(.webapp_chi(tab, correct = FALSE)),
    yates = as.numeric(.webapp_chi(tab, correct = TRUE)),
    fisher = as.numeric(pf),
    likelihood_ratio = as.numeric(.webapp_lr(tab)),
    n1_chisq = as.numeric(.webapp_cmh(tab)),
    fisher_midp = as.numeric(.webapp_midp_exact(tab)),
    fisher_midp_sas = as.numeric(.webapp_midp_classic(tab, pf))
  )
}

# The seven tests the web app implements (recalc-equivalent names).
reappraised_webapp_methods <- c(
  "pearson",
  "yates",
  "fisher",
  "likelihood_ratio",
  "n1_chisq",
  "fisher_midp",
  "fisher_midp_sas"
)
