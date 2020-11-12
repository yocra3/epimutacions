#' Function to annotate DMR resulting from epimutacions package
#' 
#' This function annotates a methylated region
#' @param data DataFrame-like object.
#' @param db Database to use for annotation. E.g: IlluminaHumanMethylationEPICanno.ilm10b2.hg19.
#' @param split Separator for CpG ids. Default ','.
#' @param epi_col CpG ids, should be rownames in database. 
#' @param gene_col Column name from where to extract gene names. Default: 'GencodeBasicV12_NAME'.
#' @param feat_col Column name from where to extract CpG feature groups. Default:  'Regulatory_Feature_Group'.
#' @param relat_col Column name from where to extract relation to island info. Default:  'Relation_to_Island'.
#' @param build Build for bioMart. Default '37'
#' @param omim Bool, if TRUE will annotate OMIMs as well. Takes a bit longer. Default TRUE
#' @return DataFrame-like object annotated.
#' @export
function(data, db, split=',', 
		# illumina annotation parameters
		epi_col='CpG_ids', gene_col='GencodeBasicV12_NAME', feat_col='Regulatory_Feature_Group', relat_col='Relation_to_Island',
		# biomart parameters
		build="37", omim=T)
{
	
	
	# get illumina data base
	message('Loading DB')
	anno <- minfi::getAnnotation(db)
	
	epids_list <- data[[epi_col]]
	
	epids_list <- strsplit(x = epids_list, split=split)
	
	
	# Basic annotation
	message('Annotating:')
	
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
	
	# OMIM annotation
	if (omim==T){
	# Using biomart  to extract OMIMs
	message('- Extracting and annotating OMIMs')
	gene_mart <- biomaRt::useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", 
									 GRCh = build, 
									 host = "www.ensembl.org", 
									 dataset="hsapiens_gene_ensembl")
	
	annotated_omim <- lapply(annotated_genes, function(gene_list){
		gene_list <- unlist(strsplit(x = gene_list, split='///'))
		# extract OMIM info from db
		omims <- biomaRt::getBM(mart = gene_mart,
								attributes = c("hgnc_symbol", "mim_morbid_accession", "mim_morbid_description"),
								filters = "hgnc_symbol",
								values = gene_list
								)
		# get accesion per gene
		accesions <- lapply(gene_list, function(gene){
			acc <- omims[omims$hgnc_symbol==gene,]$mim_morbid_accession
			acc[is.na(acc)] <- ''
			paste(acc, collapse = '/')
		})
		
		# get description per gene
		descriptions <- lapply(gene_list, function(gene){
			descr <- omims[omims$hgnc_symbol==gene,]$mim_morbid_accession
			descr[is.na(descr)] <- ''
			paste(descr, collapse = '/')
		})
		
		# per gene set
		omims_acc <- paste(accesions, collapse = '///')
		omims_desc <- paste(descriptions, collapse='///')
		return(list('OMIM_ACC'=omims_acc, 'OMIM_DESC'=omims_desc))
		
	})
	# Appending OMIM data
	annotated_omim <- do.call('rbind', annotated_omim)
	data <- cbind(data, annotated_omim)
	}
	
	return(data)
}

