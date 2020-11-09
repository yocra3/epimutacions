#' @export
EpiMutations<-function(cases, controls, num.cpgs = 10, pValue.cutoff = 0.01, 
                       cutoff =0.1, outlier.score = 0.5, 
                       nsamp = "deterministic",method = "manova")
{
  
  #Correct parameter verification
  
  if(is.null(cases))
  {
    stop("'Diseases' parameter must be introduced")
  }
  
  #Diseases length(diseases) ==1
  num_cases <- ncol(cases)
  num_controls <- ncol(controls)
  
  if(num_cases != 1) {
    stop('Number of cases has to be one')
  }
  
  if(num_controls < 2) {
    stop('Number of samples in controls (aka.reference panel) must be greater than 2')
  }
  
  if(length(method) !=1) {
    stop('Only one "method" can be chosen at a time')
  }
  
  #dataset type (GenomicRatioSet or ExpressionSet)
  type <- charmatch(class(cases), c("GenomicRatioSet", "ExpressionSet"))
  
  if(is.na(type)) {
    stop("The data type must be 'GenomicRatioSet' or 'ExpressionSet'")  
  }
  
  #Select one of the available methods
  
  selected_method <- charmatch(method, c("manova", "mlm", "iso.forest", "Mahdist.MCD"))
  
  if(is.na(selected_method)) {
    stop("The selected method must be 'manova', 'mlm','iso.forest','Mahdist.MCD'")  
  }
  
  
  #Combine control panel with the disease sample
  if(type == 1) { # GenomicRatioSets
    set <- minfi::combineArrays(controls, cases,
                                outType = c("IlluminaHumanMethylation450k",
                                            "IlluminaHumanMethylationEPIC",
                                            "IlluminaHumanMethylation27k"),
                                verbose = TRUE)
  } else if (type == 2) { #ExpressionSet
    set <- a4Base::combineTwoExpressionSet(es.control.panel,
                                           diseases)
    exprs.mat<-Biobase::exprs(set)
    fdata<-Biobase::fData(set)
  }
  
  #Obtain Phenotypic data  
  pdata <- Biobase::pData(set)
  sample <- colnames(cases)
  
  #create a variable 0,0,0,0...0,0,1  
  pdata$samp <- pdata$sampleID == sample
  
  #Create the model matrix
  model <- stats::model.matrix(~samp, pdata)
  
  #Bumphunter function from bumphunter package
  #GenomicRatioSet    
  if(type == 1) {
    bumps <- bumphunter::bumphunter(set, model, cutoff = 0.1)$table
  } else if (type == 2) { #ExpressionSet 
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
    for(i in 1:nrow(bumps)) {
      #Find beta value matrix for each bump
      beta.values <- get_betas(bumps[i, ], set)
      #manova
      if(method == "manova")
      {
        if(ncol(beta.values) > nrow(beta.values) - 2 ) {
          stop("Not enoght samples to run MANOVA")
        } else {
          bumps$manova[i]<-EpiMANOVA(beta.values, model)
        }
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