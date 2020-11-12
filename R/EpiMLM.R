#'
#'  @export
epiMLM<-function(beta.values, model)
{
  #select "model" variable columns (unique(model[,i])!= 1)
  mod <- mlm::mlm(beta.values ~ model[,2])
  statistics<-mod$aov.tab[1, c("Pr(>F)", "R2", "Pr(>F)")]
  
  return(statistics)
}