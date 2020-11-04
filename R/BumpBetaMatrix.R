#' @export

BumpBetaMatrix<-function(bump,set)
{
  #dataset type
  
  type <- charmatch(class(set), c("GenomicRatioSet", "ExpressionSet"))
  
  if(is.na(type))
    stop("The data type must be 'GenomicRatioSet' or 'ExpressionSet'")
  
  #obtain beta values matrix
  
  if(type == 1)
  {
    beta.values <- minfi::getBeta(set[bump$indexStart:bump$indexEnd, ]) 
  }
  
  else if (type == 2)
  {
    beta.values <- Biobase::exprs(set[bump$indexStart:bump$indexEnd, ])
  }
  
  beta.values<-t(na.omit(beta.values))
  
  return(beta.values)
}
