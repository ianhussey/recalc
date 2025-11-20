fisher_test_from_counts <- function(a, b, c, d, alternative = "two.sided") {
  # Build the 2x2 contingency table
  tab <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
                dimnames = list(Group = c("Row1", "Row2"),
                                Outcome = c("Col1", "Col2")))
  
  # Perform Fisher's exact test
  res <- fisher.test(tab, alternative = alternative)
  return(res)
}
Example:
  r
Copy code
# Suppose:
# Group1: 10 yes, 5 no
# Group2: 3 yes, 12 no
fisher_test_from_counts(a = 10, b = 5,
                        c = 3,  d = 12)

fisher_test_from_counts <- function(a, b, c, d, alternative = "two.sided") {
  # Build the 2x2 contingency table
  tab <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
                dimnames = list(Group = c("Row1", "Row2"),
                                Outcome = c("Col1", "Col2")))
  
  # Perform Fisher's exact test
  res <- fisher.test(tab, alternative = alternative)
  return(res)
}

# e.g.
fisher_test_from_counts(a = 29, b = 11,
                        c = 27,  d = 13)