#' Identify methylated outlier regions with a Multivariate linear model
#'
#' @param betas (matrix) Beta values of shape (n_cpgs, n_samples).
#' @param Model Formula describing the model to be fitted.
#'
#' @return A named numeric of summary statistics returned by \code{\link[mlm]{mlm}}.
#' \describe{
#' \item{F value}{The F-statistic.}
#' \item{R2}{The partial R-squared statistic.}
#' \item{Pr(>F)}{The P value.}
#' }
#' @export
epiMLM<-function(betas, model)
{
  mod <- mlm::mlm(betas ~ model[,2])
  output <- mod$aov.tab[1, c("F value", "R2", "Pr(>F)")]
  
  return(output)
}