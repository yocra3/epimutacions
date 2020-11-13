#' Return epimutations as per Aref-Eshghi et al, 2020.
#'
#' @param cases (GenomicRatioSet, ExpressionSet) Case dataset.
#' @param controls (GenomicRatioSet, ExpressionSet) Control dataset. Optional.
#' @param sample_ids (character vector) The column names in cases to compute 
#' epimutations for.
#' If missing, computes for all samples in cases.
#' @param cases_as_controls (bool) If True, all remaining cases are added in controls.
#' If False, they are ignored.
#' @param stats (string) The epimutation statistical scoring method. Choose from 
#' "manova", "mlm", "iso.forest", "Mahdist.MCD".
#' @param epis_kept (string) The epimutation filtering criteria. Either
#' "all", "aref-eshghi" or "barbosa". See \code{\link{filter_epis}}.
#' @param reduced_output (bool) If True (default), return the columns specified
#' in \code{Value}. If False, add also the columns returned by
#' \link[bumphunter]{bumphunter}.
#' @param args.bumphunter (list) Additional arguments to pass to 
#' \link[bumphunter]{bumphunter}.
#' @param args.stats (list) Additional arguments to pass to the statistical
#' scoring method called by the \code{stats} argument.
#'
#' @return A tibble of epimutation regions for sample_ids.
#' \describe{
#' \item{sample}{(string) Sample_id in colnames(cases)}
#' \item{chr}{(string) Chromosome}
#' \item{start}{(integer) Start position on chromosome}
#' \item{end}{(integer) End position on chromosome}
#' \item{length}{(integer) Region length on chromosome}
#' \item{n_cpgs}{(integer) Number of CpGs spanned by region}
#' \item{cpg_ids}{(string) Comma-separated CpG_ids in rownames(cases) spanned by region}
#' \item{beta_diff}{(numeric) Average regional methylation difference between
#' case and controls.}
#' \item{outlier_method}{(string) The outlier scoring method. Either
#' "manova", "mlm", "iso.forest" or "Mahdist.MCD".}
#' \item{outlier_score}{(numeric) The outlier score. A F-statistic for methods 
#' "manova" or "mlm". A score produced by \code{\link[isotree]{isolation.forest}} 
#' for method "iso.forest". A boolean for method "Mahdist.MCD".}
#' \item{outlier_significance}{(numeric) The outlier significance statistic.
#' P(>F-statistic) for methods "manova" or "mlm". NA for all other methods.}
#' }
#' @export
#'
#' @examples
#' data("genomicratioset") # load toy dataset
#' epi <- epimutations(genomicratioset)
epimutations <- function(
  cases,
  controls,
  sample_ids,
  cases_as_controls = T,
  stats = "manova",
  epis_kept = "all",
  reduced_output = T,
  args.bumphunter = list(cutoff=0.1),
  args.stats = list()
) {
  
  check_params(cases, controls, stats, cases_as_controls, sample_ids)

  set <- set_concat(cases, controls)
  
  if(missing(sample_ids)) sample_ids <- colnames(cases)
  epis <- lapply(
    sample_ids,
    function(sample_id) {
      epimutations_per_sample(
        set, sample_id, cases_as_controls, stats, reduced_output,
        args.bumphunter, args.stats
      )
    }
  )
  epis <- do.call(rbind, epis)
  
  epis <- filter_epis(epis, epis_kept)
  
  return(epis)
}

#' Filter epimutations.
#'
#' @param epis (data.frame) A dataframe of epimutations.
#' @param epis_kept (string) The epimutation filtering criteria. Either
#' "all", "aref-eshghi" or "barbosa".
#'
#' @return A dataframe of epimutations.
#' @export
#'
#' @examples
filter_epis <- function(epis, epis_kept){
	if(epis_kept == "aref-eshghi"){
		epis <- subset(epis,
			n_cpgs >= 3 &
				outlier_score >= 20 &
				( outlier_score >= 40 | abs(beta_diff) >= 0.2)
		)
	} else if(epis_kept == "barbosa"){
		### implement barbosa filtering criteria
	}
	return(epis)
}

epimutations_per_sample <- function(
  set,
  sample_id,
  cases_as_controls = T,
  stats = "manova",
  reduced_output = T,
  args.bumphunter = list(cutoff=0.1),
  args.stats = list()
){
  set <- filter_set(set, sample_id, cases_as_controls)
  
  design <- make_bumphunter_design(set, sample_id)
  bumps <- do.call(run_bumphunter,
                   c(list(set=set, design=design), args.bumphunter))
  
  bumps <- compute_bump_outlier_scores(set, bumps, stats, sample_id, design, args.stats)
  
  epi <- format_bumps(bumps, set, sample_id, stats, reduced_output)
  
  return(epi)
}


check_params <- function(cases, controls, stats, cases_as_controls, sample_ids){

  if(missing(cases)) {
    stop("Cases argument is missing.")
  }
  if(missing(controls) && !cases_as_controls){
    stop("If controls is missing, cases_as_controls must be set to TRUE.")
  }
  if(!class(cases) %in% c("GenomicRatioSet", "ExpressionSet")) {
    stop("Cases must be of class 'GenomicRatioSet' or 'ExpressionSet'")  
  }
  if(missing(controls) && ncol(cases) < 3){
    stop("If controls is missing, cases must contain at least 3 samples.")
  }
  if(!missing(controls) && ncol(cases) == 0){
    stop("Cases must contain at least 1 sample.")
  }
  if(!missing(controls) && !class(controls) %in% c("GenomicRatioSet", "ExpressionSet")) {
    stop("Controls must be of class 'GenomicRatioSet' or 'ExpressionSet'")  
  }
  if(!missing(controls) && ncol(controls) < 2) {
    stop("Controls must contain at least 2 samples.")
  }
  if(length(stats) != 1) {
    stop("Only one statistical method can be chosen at a time.")
  }
  if(!stats %in% c("manova", "mlm", "iso.forest", "Mahdist.MCD")) {
    stop("Statistical method must be 'manova', 'mlm','iso.forest','Mahdist.MCD'")  
  }
  if(!missing(sample_ids) && any(!sample_ids %in% colnames(cases))){
    stop("Sample_ids must match column names in the cases object.")  
  }
    n_controls <- ifelse(missing(controls), 0, ncol(controls))
    if(cases_as_controls) n_controls <- n_controls + ncol(cases) - 1
    if(stats == "manova" && n_controls < 10){
        warning("At least 10 control samples are recommended for manova statistics.")
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
#' @return The object returned by \link[bumphunter]{bumphunter}, or if no bumps
#' are found, a mock bumps object with zero rows returned by \code{empty_bumps}.
#'
#' @examples
run_bumphunter <- function(set, design, ...){
	if(class(set) == "GenomicRatioSet") {
		bumps <- bumphunter::bumphunter(set, design, ...)$table
	} else if (class(set) == "ExpressionSet") {
		bumps <- bumphunter::bumphunter(
		    object = Biobase::exprs(set),
		    design = design,
		    pos = Biobase::fData(set)$RANGE_START,
		    chr = Biobase::fData(set)$CHR,
		    ...
		)$table
	}
	if(length(bumps) == 1 && is.na(bumps)){ # bumphunter found no bumps
		bumps <- empty_bumps()
	}
	return(bumps)
}

filter_bumps <- function(bumps, min_cpgs_per_bump){
  bumps <- subset(bumps, L >= min_cpgs_per_bump)
  return(bumps)
}

compute_bump_outlier_scores <- function(set, bumps, stats, sample, model, args.stats){
  bumps$outlier_score <- bumps$outlier_significance <- rep(NA_real_, nrow(bumps))
  bumps$beta_diff <- rep(NA_real_, nrow(bumps))
  for(i in seq_len(nrow(bumps))) {
    betas <- get_betas(bumps[i, ], set)
    if(stats == "manova") {
        manova_has_enough_samples <- nrow(betas) >= ncol(betas) + 2
      if(manova_has_enough_samples) {
          stats_manova <- epi_manova(betas, model, sample)
          bumps$outlier_score[i] <- stats_manova["approx F"]
          bumps$outlier_significance[i] <- stats_manova["Pr(>F)"]
          bumps$beta_diff[i] <- stats_manova["beta_mean_abs_diff"]
      }
    } else if(stats == "mlm") {
        stats_mlm <- epiMLM(betas, model)
      bumps$outlier_score[i] <- stats_mlm["F value"]
      bumps$outlier_significance[i] <- stats_mlm["Pr(>F)"]
    } else if(stats == "iso.forest") {
      bumps$outlier_score[i] <- epiIsolationForest(betas, sample)
    } else if(stats == "Mahdist.MCD") {
    	bumps$outlier_score[i] <- do.call(
    		epiMahdist.MCD, c(list(betas=betas, sample=sample), args.stats))
    }
  }
  
  return(bumps)
}

select_outlier_bumps <- function(bumps, stats, pValue.cutoff, fStat_min, betaDiff_min, outlier.score){

    if(stats == "manova"){
        outliers <- subset(
            bumps,
            outlier_score >= fStat_min &
                beta_diff >= betaDiff_min &
                outlier_significance < pValue.cutoff 
        )
    } else if(stats == "mlm"){
        outliers <- subset(bumps, outlier_significance < pValue.cutoff)
    } else if(stats == "iso.forest"){
        outliers <- subset(bumps, outlier_score > outlier.score)
    } else if(stats == "Mahdist.MCD"){
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
format_bumps <- function(bumps, set, sample, stats, reduced){
	df_out <- tibble::as_tibble(bumps)
	df_out$sample <- df_out$outlier_method <- character(nrow(df_out))
	if(nrow(bumps) > 0){
		df_out$sample <- sample
		df_out$outlier_method <- stats
	}
	df_out$length <- df_out$end - df_out$start + 1
	df_out$n_cpgs <- df_out$indexEnd - df_out$indexStart + 1
	df_out$cpg_ids <- mapply(
		function(rown, i_st, i_end){paste(rown[i_st:i_end], collapse=",")},
		df_out$indexStart,
		df_out$indexEnd,
		MoreArgs = list(rown = rownames(set))
	)
	if(reduced){
		reduced_col <- c("sample", "chr", "start", "end", "length", "n_cpgs", "cpg_ids",
		                 "beta_diff", "outlier_method", "outlier_score", "outlier_significance")
		df_out <- df_out[, reduced_col]
	}
	return(df_out)
}

empty_bumps <- function(){
	bumps <- tibble::tibble(
		"chr" = character(),
		"start" = integer(),
		"end" = integer(),
		"value" = numeric(),
		"area" = numeric(),
		"cluster" = numeric(),
		"indexStart" = integer(),
		"indexEnd" = integer(),
		"L" = numeric()
	)
	return(bumps)
}