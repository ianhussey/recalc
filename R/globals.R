#' @importFrom stats setNames
NULL

# Column names used in non-standard evaluation (dplyr verbs, ggplot2 aes())
# are not visible to R CMD check; declare them here to silence the note.
utils::globalVariables(c(
  ".data",
  "min_d",
  "max_d",
  "min_d_rounded",
  "max_d_rounded",
  "min_p",
  "max_p",
  "min_p_rounded",
  "max_p_rounded",
  "d_inbounds",
  "d_inbounds_hull",
  "p_inbounds",
  "p_inbounds_hull"
))
