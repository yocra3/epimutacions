#' @export
EpiMutations <- function(
  cases,
  controls,
  num.cpgs = 10,
  pValue.cutoff = 0.01,
  cutoff = 0.1,
  outlier.score = 0.5,
  nsamp = "deterministic",
  method = c("manova", "mlm", "iso.forest", "Mahdist.MCD")
) {
  
  check_params(cases, controls, method)

  set <- set_concat(cases, controls)
  
  model <- make_bumphunter_model(set, sample_id=colnames(cases))
  bumps <- run_bumphunter(set, model, cutoff)
  
  check_bumps(bumps)
  bumps <- filter_bumps(bumps, min_cpgs_per_bump=num.cpgs)
  
  bumps <- compute_bump_outlier_scores(set, bumps, method, colnames(cases), model, nsamp)
  bumps <- select_outlier_bumps(bumps, method, pValue.cutoff, outlier.score)
  
  return(format_bumps(bumps, colnames(cases)))
}

check_params <- function(cases, controls, method){
  if(is.null(cases)) {
    stop("'Diseases' parameter must be introduced")
  }
  if(ncol(cases) != 1) {
    stop('Number of cases has to be one')
  }
  if(ncol(controls) < 2) {
    stop('Number of samples in controls (aka.reference panel) must be greater than 2')
  }
  if(length(method) !=1) {
    stop('Only one "method" can be chosen at a time')
  }
  type <- charmatch(class(cases), c("GenomicRatioSet", "ExpressionSet"))
  if(is.na(type)) {
    stop("The data type must be 'GenomicRatioSet' or 'ExpressionSet'")  
  }
  selected_method <- charmatch(method, c("manova", "mlm", "iso.forest", "Mahdist.MCD"))
  if(is.na(selected_method)) {
    stop("The selected method must be 'manova', 'mlm','iso.forest','Mahdist.MCD'")  
  }
}

set_concat <- function(cases, controls){
  if(class(cases) == "GenomicRatioSet") {
    return(minfi::combineArrays(
      controls, cases,
      outType = c("IlluminaHumanMethylation450k",
        "IlluminaHumanMethylationEPIC",
        "IlluminaHumanMethylation27k"),
        verbose = TRUE)
      )
  } else if (class(cases) == "ExpressionSet") {
    return(a4Base::combineTwoExpressionSet(controls, cases))
  }
}

make_bumphunter_model <- function(set, sample_id){
  # model over samples is 0,0,0,0...0,0,1
  Biobase::pData(set)$samp <- Biobase::pData(set)$sampleID == sample_id
  return(stats::model.matrix(~samp, Biobase::pData(set)))
}

run_bumphunter <- function(set, model, cutoff){
  if(class(set) == "GenomicRatioSet") {
    return(bumphunter::bumphunter(set, model, cutoff = cutoff)$table)
  } else if (class(set) == "ExpressionSet") {
    return(
      bumphunter::bumphunter(
        object = Biobase::exprs(set),
        design = model,
        pos = Biobase::fData(set)$RANGE_START,
        chr = Biobase::fData(set)$CHR,
        cutoff = cutoff)
      $table
    )
  }
}

check_bumps <- function(bumps){
  if(is.na(bumps[1,1])) {
    stop("Bumhunter has returned NAs.")
  }
}

filter_bumps <- function(bumps, min_cpgs_per_bump){
  bumps <- subset(bumps, L >= min_cpgs_per_bump)
  return(bumps)
}

compute_bump_outlier_scores <- function(set, bumps, method, sample, model, nsamp){
  
  for(i in seq_len(nrow(bumps))) {
    beta.values <- get_betas(bumps[i, ], set)
    if(method == "manova") {
      if(ncol(beta.values) > nrow(beta.values) - 2 ) {
        stop("Not enough samples to run MANOVA")
      }
      bumps$manova[i]<-EpiMANOVA(beta.values, model)
    }
    if(method == "mlm") {
      bumps$mlm[i]<-epiMLM(beta.values, model)  
    }
    if(method == "iso.forest") {
      bumps$iso[i]<-epiIsolationForest(beta.values, sample)
    }
    if(method == "Mahdist.MCD") {
      bumps$MahMCD[i]<-epiMahdist.MCD(beta.values, nsamp, sample)
    }
  }
  return(bumps)
}

select_outlier_bumps <- function(bumps, method, pValue.cutoff, outlier.score){

  if(method == "manova"){
    outliers <- subset(bumps, manova < pValue.cutoff)
  }
  if(method == "mlm"){
    outliers <- subset(bumps, mlm < pValue.cutoff)
  }
  if(method == "iso.forest"){
    outliers <- subset(bumps, iso > outlier.score)
  }
  if(method == "Mahdist.MCD"){
    outliers <- subset(bumps, MahMCD == TRUE)
  }
  
  return(outliers)
}

format_bumps <- function(bumps, sample){
  bumps$sample <- sample
  return(tibble::as_tibble(bumps))
}