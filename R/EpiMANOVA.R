#' Identifies methylated outlier regions using manova
#' 
#' This function identifies regions with CpGs being outliers 
#' 
#' @param betas Beta values matrix
#' @param Model Formula describing the model to be fitted
#' @return F statistic
#' @export
#' 
epi_manova <-  function(betas, model){
  
  mod<-manova(betas ~ model)
  mod_summary <- summary(mod)$stats
  F_stat <- mod_summary[1,"approx F"]
  return(F_stat)
}
