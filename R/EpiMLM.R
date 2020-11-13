#' Identifies methylated outlier regions using Multivariate linear model
#' 
#' This function identifies regions with CpGs being outliers 
#' 
#' @param betas Beta values matrix
#' @param Model Formula describing the model to be fitted
#' @return F statistic, R2 test statistic and Pillai
#' @export
#' 
#' @export
epiMLM<-function(betas, model)
{
  #select "model" variable columns (unique(model[,i])!= 1)
  
 
  mod <- mlm::mlm(betas ~ model[,2])
  output<-mod$aov.tab[1, c("approx F", "R2", "Pr(>F)")]
  
  return(output)
}