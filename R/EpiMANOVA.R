#' Identifies methylated outlier regions using manova
#' 
#' This function identifies regions with CpGs being outliers 
#' 
#' @param betas Beta values matrix
#' @param Model Formula describing the model to be fitted
#' @return F statistic
#' @export
#' 
epi_manova <-  function(betas, model, sample){
  
  # Fit the manova model
  mod <- manova(betas ~ model)
  # Model summary
  mod_summary <- summary(mod)$stats
  # Obtain F statistic
  F_stat <- mod_summary[1,"approx F"]
  
  # Filter out epimutation (F-statistic > 40 or F-statistic < 20 and mean difference > 0.2)
  # Calculate the mean difference
  keep_case <- which(sample %in% rownames(betas))
  case <- betas[keep_case,]
  controls <- betas[-keep_case,]
  coltrols_mean <- colMeans(controls)
  mean_difference <- coltrols_mean-case
  check_mean_difference<-mean_difference > 0.2
  check_mean_difference<-unique(check_mean_difference)
  length_check_mean_difference<-length(check_mean_difference)
  
  if(length_check_mean_difference == 1 & isTRUE(check_mean_difference)){# F-statistic < 20
    output <- F_stat < 20
  }else{# F-statistic > 40
    output <- F_stat > 40
  }
  return(output)
}