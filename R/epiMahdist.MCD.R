#' @export
epiMahdist.MCD<-function(beta.values, nsamp = c("best", "exact", "deterministic"), sample)
{
  if(is.null(beta.values))
  {
    stop("'beta.values' data frame must be introduced")
  }
  
  if(nsamp != "best" & nsamp != "exact" & nsamp != "deterministic")
  {
    stop("'nsamp' must be 'best','exact' or 'deterministic'")
    
  }
  
  if(is.null(sample))
  {
    stop("'sample' corresponding to rare disease must be introduced")
  }
  
  if(length(sample)!= 1)
  {
    stop("Only one 'sample' with suspected rare disease can be introduced") 
  }

  
  #run the  mcd model
  mcd<-robustbase::covMcd(beta.values, nsamp = nsamp)
  #get mcd stimate of location
  mean_mcd <- mcd$center
  #get mcd estimate scatter
  cov_mcd <- mcd$cov
  #get inverse of scatter
  cov_mcd_inv <- solve(cov_mcd)
  
  #compute the robust distance
  robust_dist <- apply(beta.values, 1, function(x){
    x <- (x - mean_mcd)
    dist <- sqrt((t(x)  %*% cov_mcd_inv %*% x))
    return(dist)
  })
  
  #set cutoff using chi square distribution
  threshold <- sqrt(qchisq(p = 0.975, df = ncol(beta.values)))
  
  #find outliers
  outliers <-  which(robust_dist >= threshold)
  
  #suspected rare disease sample within the outliers
  
  rd.sample.out<-sample %in% names(outliers)
  
  return(rd.sample.out)
}