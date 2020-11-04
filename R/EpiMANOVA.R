EpiMANOVA <-  function(beta.values, model)
{
  modSum <- summary(manova(beta.values ~ model[,2]))
  
  p.value<-modSum$stats[1,"Pr(>F)"]
  
  return(p.value)
}