#' Function to annotate
#' 
#' This function annotates a methylated region
#' TxDbannotation package is required
#' @param data DataFrame-like object.
#' @param db Database to use for annotation. E.g: IlluminaHumanMethylationEPICanno.ilm10b2.hg19.
#' @param epi_col CpG ids, should be rownames in database.
#' @param gene_col Column name from where to extract gene names. Default: 'GencodeBasicV12_NAME'.
#' @param feat_col Column name from where to extract CpG feature groups. Default:  'Regulatory_Feature_Group'.
#' @param relat_col Column name from where to extract relation to island info. Default:  'Relation_to_Island'.
#' @return DataFrame-like object annotated.
#' @export
annotate_CpG <- function(data, db, epi_col='CpG_ids', gene_col='GencodeBasicV12_NAME', feat_col='Regulatory_Feature_Group', relat_col='Relation_to_Island'){
	message('Annotating:')
	anno <- getAnnotation(db)
	epids_list <- data[[epi_col]]
	message(paste('-', gene_col))
	annotated_genes <- lapply(epids_list, function(x){
		paste(anno[x,][[gene_col]], collapse='///')
	})
	message(paste('-', feat_col))
	annotated_regions <- lapply(epids_list, function(x){
		paste(anno[x,][[feat_col]], collapse='///')
	})
	message(paste('-', relat_col))
	annotated_relation <- lapply(epids_list, function(x){
		paste(anno[x,][[relat_col]], collapse='///')
	})
	
	
	data[[gene_col]] <- annotated_genes
	data[[feat_col]] <- annotated_regions
	data[[relat_col]] <- annotated_relation
	return(data)
}