
# source of equation
# https://matthewbjane.quarto.pub/pre-post-correlations/#scenario-3-is-the-t-statistic-from-a-paired-t-test-available

prepost_cor_from_dependent_ttest <- function(t, SD1, SD2, n, m_diff){
  r <- 
    (t^2 * (SD1^2 + SD2^2)) - (n * m_diff) / 
    (2 * t^2 * SD1 * SD2)
  
  return(r)
}
