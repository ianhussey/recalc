#' recalc: Recalculate Statistical Results from Reported Summary and Test Statistics
#'
#' Forensic meta-science tools that recalculate statistical results from the
#' summary and test statistics reported in a paper. Because reported numbers are
#' rounded and the analysis that produced them involves choices that are often
#' left unstated, a single reported result is compatible with a range of
#' recalculated values rather than a single one. The `recalc_*()` functions
#' recompute p-values and standardised effect sizes across the multiverse of
#' defensible choices (input/output rounding, Student's versus Welch's t-tests,
#' one- versus two-sided tests, SD/SE confusion, and so on) and report whether
#' the reported value falls inside the resulting interval.
#'
#' @seealso Useful entry points:
#'   \code{\link{recalc_independent_t}} (independent-samples t-test),
#'   \code{\link{recalc_chisq}} (contingency tables),
#'   \code{\link{recalc_regression_from_cor_mat}} (regression from a correlation
#'   matrix), \code{\link{recalc_mediation}} (mediation), and
#'   \code{\link{recalc_missing_subgroup}} (whole-vs-subgroup checks).
#'
#' @keywords internal
"_PACKAGE"
