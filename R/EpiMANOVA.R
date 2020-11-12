#' Identifies methylated outlier regions using manova
#' 
#' This function identifies regions with CpGs being outliers 
#' 
#' @param betas Beta values matrix
#' @param Model Formula describing the model to be fitted
#' @param Sample_id Character vector specifying the name of the sample
#' to compute epimutations
#' @return F statistic, multivariate test statistic and Pillai
#' @export
#' 
epi_manova <-  function(betas, model, sample_id){
  
  # Fit the manova model
  mod <- manova(betas ~ model)
  
  # Model summary
  mod_summary <- summary(mod)$stats
  
  # Obtain statistics (F statistic, pillai, p value)
  statistics<- mod_summary[1,c("approx F", "Pillai","Pr(>F)")]
  
  # Calculate the beta mean difference
  case_row <- which(rownames(betas) %in% sample_id)
  beta_mean_difference <- abs(mean(colMeans(betas[-case_row,]) - betas[case_row,]))

  output<-c(statistics, beta_mean_difference)

  return(output)
}