#' Annotate epimutations results
#' 
#' Add information about close genes and regulatory elements for epimutations.
#' 
#' @param epi_results Output from `epimutations`
#' @param db String with the Illumina annotation package used to annotate the CpGs
#' @param build Genomic build where the epimutations are mapped. By default, it is 
#' GRCh37. Set build to `NULL` to use GRCh38.
#' @param ... Further arguments passed to `annotate_cpg`
#' 
#' @return Input object `epi_results` with additional columns with information about
#' the genes or overlapping regulatory features. See `annotate_cpg` and `add_ensemble_regulatory` 
#' for an in-depth description of these variables.
#' 
#' @export
annotate_epimutations <- function(epi_results, db = "IlluminaHumanMethylationEPICanno.ilm10b2.hg19", build = "37", ...){
	
	## Add gene mapping and CpG island information
	epi_results <- annotate_cpg(epi_results, db = db,  ...)
	
	## Add ENSEMBL regulatory regions
	epi_results <- add_ensemble_regulatory(epi_results, build = build)
	epi_results
}