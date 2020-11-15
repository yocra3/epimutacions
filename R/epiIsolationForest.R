#' Identifies methylated outlier regions using isolation forest
#' 
#' This function identifies regions with CpGs being outliers 
#' 
#' @param betas (matrix) Beta values of shape (n_cpgs, n_samples).
#' @param Sample_id Character vector specifying the name of the sample
#' to compute epimutations
#' @return The outlier score for the given case sample
#' @export
#' 
#'  @export
epiIsolationForest <- function(betas, sample_id)
{

  
  # Generate train and test data frame
  betas <- as.data.frame(betas)
  train <- betas[!(row.names(betas) %in% sample_id),]
  # Check if train dataset has enought variables
  if(nrow(train) < 5){
    stop("Not enough samples to run isolation forest")
  }
  test <- betas[sample_id,]
  
  # Run the isolation forest methods 
      iso <- isotree::isolation.forest(df = train)
      # Predict
      outlier.score <- predict(iso,test)

  return(outlier.score)
}