#' @export
epiMLM<-function(beta.values, model)
{
  #select "model" variable columns (unique(model[,i])!= 1)
  
 
  mod <- mlm::mlm(beta.values ~ model[,2])
  output<-mod$aov.tab[1, c("F value", "R2", "Pr(>F)")]
  
  return(output)
}