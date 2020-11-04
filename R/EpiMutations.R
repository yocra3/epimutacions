#' @export
EpiMutations<-function(diseases, num.cpgs = 10, pValue.cutoff = 0.01, 
                       cutoff =0.1, outlier.score = 0.5, 
                       nsamp = "deterministic",method = "manova")
{
  
  #Correct parameter verification
  
  if(is.null(diseases))
  {
    stop("'Diseases' parameter must be introduced")
  }
  
  #Diseases length(diseases) ==1
  num.sample.diseases<-dim(diseases)[2]
  
  if (num.sample.diseases !=1)
  {
    stop("'diseases' parameter sample number must be 1")
  }
  
  if (length(method)!=1)
  {
    stop(" Only one 'method' can be chosen at a time")
  }
  
  #dataset type (GenomicRatioSet or ExpressionSet)
  type <- charmatch(class(diseases), c("GenomicRatioSet", "ExpressionSet"))
  
  if(is.na(type))
  {
    stop("The data type must be 'GenomicRatioSet' or 'ExpressionSet'")  
  }
  
  #Select one of the available methods
  
  method.selected<-charmatch(method,c("manova","mlm","iso.forest","Mahdist.MCD"))
  
  if(is.na(method.selected))
  {
    stop("The selected method must be 'manova', 'mlm','iso.forest','Mahdist.MCD'")  
  }
  
  
  #Combine control panel with the disease sample
  #GenomicRatioSets
  if(type == 1)
  {
    
    set <- minfi::combineArrays(grs.control.panel, diseases,
                                outType = c("IlluminaHumanMethylation450k",
                                            "IlluminaHumanMethylationEPIC",
                                            "IlluminaHumanMethylation27k"),
                                verbose = TRUE)
  }
  
  #ExpressionSet   
  else if (type == 2)
  {
    
    set <- a4Base::combineTwoExpressionSet(es.control.panel,
                                           diseases)
    
    exprs.mat<-Biobase::exprs(set)
    fdata<-Biobase::fData(set)
  }
  
  #Obtain Phenotypic data  
  pdata <- Biobase::pData(set)
  sample<-colnames(diseases)
  
  #create a variable 0,0,0,0...0,0,1  
  pdata$samp <- pdata$sampleID == sample
  
  #Create the model matrix
  model <- stats::model.matrix(~ samp, pdata)
  
  #Bumphunter function from bumphunter package
  #GenomicRatioSet    
  if(type == 1)
  {
    bumps <- bumphunter::bumphunter(set, model, cutoff = 0.1)$table
  }
  
  #ExpressionSet 
  else if (type == 2)
  {
    bumps <- bumphunter::bumphunter(object = exprs.mat,
                                    design = model,
                                    pos = fdata$RANGE_START,
                                    chr = fdata$CHR,
                                    cutoff = 0.1)$table
  }
  
  #Outlier identification using multiple methods
  
  if(!is.na(bumps[1,1]))
  {
    
    bumps$sample <- sample
    
    #delect bumps with at least selected "num.cpgs"
    bumps <- subset(bumps, L >= num.cpgs)
    
    #Find beta value matrix for each bump
    
    #Outlier identification using multiple statistical approach 
    for( i in 1:dim(bumps)[1])
    {
      #Find beta value matrix for each bump
      beta.values<-BumpBetaMatrix(bumps[i,],set)
      #manova
      if(method == "manova")
      {
        bumps$manova[i]<-EpiMANOVA(beta.values,model) 
      }
      #mlm
      if(method == "mlm")
      {
        bumps$mlm[i]<-epiMLM(beta.values,model)  
      }
      if(method == "iso.forest")
      {
        bumps$iso[i]<-epiIsolationForest(beta.values, sample)
        
      }
      if(method == "Mahdist.MCD")
      {
        bumps$MahMCD[i]<-epiMahdist.MCD(beta.values, nsamp, sample)
      }
    }
    
    #Subset bumps using p value
    
    if(method == "manova")
    {
      outliers.epi.mutations <- subset(bumps, manova < pValue.cutoff)
    }
    if(method == "mlm")
    {
      outliers.epi.mutations <- subset(bumps, mlm < pValue.cutoff)
    }
    if(method == "iso.forest")
    {
      outliers.epi.mutations <- subset(bumps, iso > outlier.score)
      
    }
    if(method == "Mahdist.MCD")
    {
      outliers.epi.mutations <- subset(bumps, MahMCD == TRUE)
    }
  }
  
  outliers.epi.mutations<-tibble::as_tibble(outliers.epi.mutations)
  return(outliers.epi.mutations)
}