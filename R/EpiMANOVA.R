#' Identifies methylated outlier regions using manova
#' 
#' This function identifies regions with CpGs being outliers 
#' 
#' @param betas (matrix) Beta values of shape (n_cpgs, n_samples).
#' @param Model Formula describing the model to be fitted.
#' @param Sample_id (string) The sample in colnames(betas)
#' to compute epimutations on.
#' @return A named numeric of summary statistics.
#' \describe{
#' \item{approx F}{The approximate F-statistic returned by \code{\link[stats]{summary.manova}}}
#' \item{Pillai}{The Pillai statistic returned by \code{\link[stats]{summary.manova}}}
#' \item{Pr(>F)}{The P value returned by \code{\link[stats]{summary.manova}}}
#' \item{beta_mean_abs_diff}{The mean absolute difference in beta values between sample and controls.}
#' }
#' @export
#' 
epi_manova <- function(betas, model){
  
  # Fit the manova model
  mod <- manova(betas ~ model)
  mod_summary <- summary(mod)$stats
  
  # Obtain statistics (F statistic, Pillai statistic, p value)
  statistics <- mod_summary[1, c("approx F", "Pillai","Pr(>F)")]

  return(statistics)
}