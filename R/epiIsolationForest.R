#' @export
epiIsolationForest<-function(beta.values, disease.sample.name)
{
  if(is.null(beta.values))
  {
   stop("'beta.values' data frame must be introduced")
  }
  
  if(is.null(disease.sample.name) | class(disease.sample.name) != "character")
  {
   stop("You must provide the name of the disease in 'disease.sample.name'")
  }
  
  #Generate train and test(sample with suspected disease) data frame
  beta.values<-as.data.frame(beta.values)
  train<-beta.values[!(row.names(beta.values) %in% disease.sample.name),]
  test<-beta.values[disease.sample.name,]
  
  #Run the isolation forest methods 
      iso<-isotree::isolation.forest(df = train)
      #Predict
      outlier.score<-predict(iso,test)

  return(outlier.score)
}