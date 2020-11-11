#' Return epimutations as per Aref-Eshghi et al, 2020.
#'
#' @param cases (GenomicRatioSet, ExpressionSet) Case dataset.
#' @param controls (GenomicRatioSet, ExpressionSet) Control dataset. Optional.
#' @param sample_ids (character vector) The column names in cases to compute 
#' epimutations for.
#' If missing, computes for all samples in cases.
#' @param cases_as_controls (bool) If True, all remaining cases are added in controls.
#' If False, they are ignored.
#' @param args.bumphunter (list) Additional arguments to pass to 
#' \link[bumphunter]{bumphunter}.
#' @param num.cpgs (integer) Bumps containing less cpgs than num.cpgs are discarded.
#' @param pValue.cutoff 
#' @param outlier.score 
#' @param nsamp 
#' @param method (string) The outlier scoring method. Choose from 
#' "manova", "mlm", "iso.forest", "Mahdist.MCD".
#'
#' @return A tibble of epimutation regions for sample_id.
#' \describe{
#' \item{sample}{(string) Sample_id in colnames(cases)}
#' \item{chr}{(string) Chromosome}
#' \item{start}{(integer) Start position on chromosome}
#' \item{end}{(integer) End position on chromosome}
#' \item{cpg_ids}{(string) Comma-separated CpG_ids in rownames(cases) spanned by region}
#' \item{outlier_method}{(string) The outlier scoring method. Either
#' "manova", "mlm", "iso.forest" or "Mahdist.MCD".}
#' \item{outlier_score}{(numeric) The outlier score. A p-value for methods 
#' "manova" or "mlm". A score produced by \code{\link[isotree]{isolation.forest}} 
#' for method "iso.forest". A boolean for method "Mahdist.MCD". }
#' }
#' @export
#'
#' @examples
#' data("genomicratioset") # load toy dataset
#' epi <- epimutations(
#'   genomicratioset,
#'   num.cpgs = 2,
#'   method = "manova"
#' )
epimutations <- function(
  cases,
  controls,
  sample_ids,
  cases_as_controls = T,
  args.bumphunter = list(cutoff=0.1),
  num.cpgs = 10,
  pValue.cutoff = 0.01,
  outlier.score = 0.5,
  nsamp = "deterministic",
  method = c("manova", "mlm", "iso.forest", "Mahdist.MCD"),
  reduced_output = T
) {
  
  check_params(cases, controls, method, cases_as_controls)

  set <- set_concat(cases, controls)
  
  if(missing(sample_ids)){
    sample_ids <- colnames(cases)
  }
  epis <- lapply(
    sample_ids,
    function(sample_id) {
      epimutations_per_sample(
        set, sample_id, cases_as_controls, args.bumphunter, num.cpgs,
        pValue.cutoff, outlier.score, nsamp, method, reduced_output
      )
    }
  )
  epis <- do.call(rbind, epis)
  
  return(epis)
}

epimutations_per_sample <- function(
  set,
  sample_id,
  cases_as_controls = T,
  args.bumphunter = list(cutoff=0.1),
  num.cpgs = 10,
  pValue.cutoff = 0.01,
  outlier.score = 0.5,
  nsamp = "deterministic",
  method = c("manova", "mlm", "iso.forest", "Mahdist.MCD"),
  reduced_output = T
){
  set <- filter_set(set, sample_id, cases_as_controls)
  
  design <- make_bumphunter_design(set, sample_id)
  bumps <- do.call(run_bumphunter,
                   c(list(set=set, design=design), args.bumphunter))
  
  check_bumps(bumps)
  bumps <- filter_bumps(bumps, min_cpgs_per_bump=num.cpgs)
  
  bumps <- compute_bump_outlier_scores(set, bumps, method, sample_id, design, nsamp)
  bumps <- select_outlier_bumps(bumps, method, pValue.cutoff, outlier.score)
  epi <- format_bumps(bumps, set, sample_id, method, reduced_output)
  
  return(epi)
}

check_params <- function(cases, controls, method, cases_as_controls){
  if(is.null(cases)) {
    stop("'Diseases' parameter must be introduced")
  }
  if(missing(controls) && !cases_as_controls){
    stop("If controls is missing, cases_as_controls must be set to TRUE.")
  }
  if(!missing(controls) && ncol(controls) < 2) {
    stop('Number of samples in controls (aka.reference panel) must be greater than 2')
  }
  if(length(method) != 1) {
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
  Biobase::pData(cases)[["origin"]] <- "case"
  if(missing(controls)){
    return(cases)
  }
  Biobase::pData(controls)[["origin"]] <- "control"
  if(class(cases) == "GenomicRatioSet") {
    set <- minfi::combineArrays(
      controls, cases,
      outType = c("IlluminaHumanMethylation450k",
        "IlluminaHumanMethylationEPIC",
        "IlluminaHumanMethylation27k"),
        verbose = TRUE
    )
  } else if (class(cases) == "ExpressionSet") {
    set <- a4Base::combineTwoExpressionSet(controls, cases)
  }
  return(set)
}

filter_set <- function(set, sample_id, cases_as_controls){
  if(!cases_as_controls){
    mask <- (
      Biobase::pData(set)[["origin"]] == "control" 
        | colnames(set) == sample_id
    )
    set <- set[, mask]
  }
  return(set)
}

make_bumphunter_design <- function(set, sample_id){
  # model is single sample, no covariates: 0,0,0,0...0,0,1
  Biobase::pData(set)$samp <- colnames(set) == sample_id
  return(stats::model.matrix(~samp, Biobase::pData(set)))
}

#' Run the Bumphunter algorithm.
#'
#' See \link[bumphunter]{bumphunter} doc for details.
#'
#' @keywords internal
#' 
#' @param set A GenomicRatioSet or ExpressionSet
#' @param design Design matrix with rows representing samples and columns 
#' representing covariates. Regression is applied to each row of set.
#' @param ... Extra arguments to pass to \link[bumphunter]{bumphunter} 
#'
#' @return The object returned by \link[bumphunter]{bumphunter}
#'
#' @examples
run_bumphunter <- function(set, design, ...){
  if(class(set) == "GenomicRatioSet") {
    return(bumphunter::bumphunter(set, design, ...)$table)
  } else if (class(set) == "ExpressionSet") {
    return(
      bumphunter::bumphunter(
        object = Biobase::exprs(set),
        design = design,
        pos = Biobase::fData(set)$RANGE_START,
        chr = Biobase::fData(set)$CHR,
        ...
      )$table
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
  bumps$outlier_score <- character(nrow(bumps))
  for(i in seq_len(nrow(bumps))) {
    beta.values <- get_betas(bumps[i, ], set)
    if(method == "manova") {
      if(ncol(beta.values) > nrow(beta.values) - 2 ) {
        stop("Not enough samples to run MANOVA")
      }
      bumps$outlier_score[i] <- EpiMANOVA(beta.values, model)
    } else if(method == "mlm") {
      bumps$outlier_score[i] <- epiMLM(beta.values, model)  
    } else if(method == "iso.forest") {
      bumps$outlier_score[i] <- epiIsolationForest(beta.values, sample)
    } else if(method == "Mahdist.MCD") {
      bumps$outlier_score[i] <- epiMahdist.MCD(beta.values, nsamp, sample)
    }
  }
  
  return(bumps)
}

select_outlier_bumps <- function(bumps, method, pValue.cutoff, outlier.score){

  if(method == "manova" || method == "mlm"){
    outliers <- subset(bumps, outlier_score < pValue.cutoff)
  } else if(method == "iso.forest"){
    outliers <- subset(bumps, outlier_score > outlier.score)
  } else if(method == "Mahdist.MCD"){
    outliers <- subset(bumps, outlier_score == TRUE)
  }
  
  return(outliers)
}

#' Convert bumps object to a tibble with additional formatting.
#'
#' @keywords internal
#'
#' @param bumps (bumps) Bump object returned by \link[bumphunter]{bumphunter}.
#' @param set (GenomicRatioSet, ExpressionSet) The dataset on which the bumps
#' were computed.
#' @param sample (string) The sample name for all the bumps.
#' @param reduced (bool) If True, only return a limited number of columns.
#' If False, return a full set of columns.
#'
#' @return A tibble dataframe where each row is an epimutation.
#'
#' @examples
format_bumps <- function(bumps, set, sample, method, reduced){
  df_out <- tibble::as_tibble(bumps)
  df_out$sample <- df_out$outlier_method <- character(nrow(df_out))
  if(nrow(bumps) > 0){
    df_out$sample <- sample
    df_out$outlier_method <- method
  }
  df_out$cpg_ids <- mapply(
    function(rown, i_st, i_end){paste(rown[i_st:i_end], collapse=",")},
    df_out$indexStart,
    df_out$indexEnd,
    MoreArgs = list(rown = rownames(set))
  )
  if(reduced){
    reduced_col <- c("sample", "chr", "start", "end", "cpg_ids",
                     "outlier_method", "outlier_score")
    df_out <- df_out[, reduced_col]
  }
  return(df_out)
}