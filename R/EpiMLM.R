#' @export
epiMLM<-function(beta.values, model)
{
  #select "model" variable columns (unique(model[,i])!= 1)
  
 
  mod <- mlm::mlm(beta.values ~ model[,2])
  p.value<-mod$aov.tab[1, "Pr(>F)"]
  
  return(p.value)
}