#' Add ENSEMBL regulatory regions to epimutations
#' 
#' @param epimutations Result from `EpiMutations`
#' @param build Build used to define epimutations coordinates. By default, it is 37,
#' corresponding to Illumina annotation.
#' 
#' @return Object of `epimutations` with some additional variables describing regulatory
#' elements from ENSEMBL. Note that a single epimutation might overlap with more than one
#' regulatory region. In that case, the different regulatory regions are separated by `///`.
#' \itemize{
#'  \item{ensembl_reg_id}{Region identifier from ENSEMBL}
#'  \item{ensembl_reg_coordinates}{Coordinates for the ENSEMBL regulatory regions}
#'  \item{ensembl_reg_type}{Type of regulatory region}
#'  \item{ensembl_reg_tissues}{Activity of the regulatory region per tissue. The different
#'  activation states are separated by `/`}
#' }
#'
add_ensemble_regulatory <- function(epimutations, build = "37"){
	
	## Create connection to ENSEMBL 
	mart <- biomaRt::useEnsembl(biomart = "regulation", GRCh = build)
	ensembl <- biomaRt::useDataset(dataset = "hsapiens_regulatory_feature", 
								   mart = mart)

	reg_res <- lapply(seq_len(nrow(epimutations)), function(i) {
		get_ENSEMBL_data(epimutations[i, "chromosome"], 
						 epimutations[i, "start"],
						 epimutations[i, "end"], mart = ensembl)
	})
	reg_res_df <- Reduce(rbind, reg_res)
	cbind(epimutations, reg_res_df)

}

#' Get ENSEMBL regulatory features overlapping a genomic region
#' 
#' This function queries for ENSEMBL regulatory features and collapse them to return
#' a single record.
#' 
#' @param chromosome Chromosome of the region
#' @param start Start of the region
#' @param end End of the region
#' @param mart `Mart` object to perform the ENSEMBL query
#' @return `data.frame` of one row with the ENSEMBL regulatory regions overlapping
#' the genomic coordinate
get_ENSEMBL_data <- function(chromosome, start, end, mart){
	bm <- biomaRt::getBM(attributes = c("activity", "regulatory_stable_id", 
								  "chromosome_name", "chromosome_start",
								  "chromosome_end", "feature_type_name",
								  "epigenome_name"), 
				   filters = c('chromosome_name','start','end'),
				   values = list(chromosome, start, end),
				   mart = mart)
	out_ens <- process_ENSEMBL_results(bm)
	out_ens
}

#' Process data from ENSEMBL 
#' 
#' Process data from ENSEMBL to combine results from the same regulatory elements in 
#' a unique record.
#' 
#' @param ensembl_res Results from `biomaRt::getBM`
#' @return `data.frame` of one row after collapsing the input ENSEMBL regulatory regions 
process_ENSEMBL_results <- function(ensembl_res){
	
	reg_elements <- unique(ensembl_res$regulatory_stable_id)
	reg_vals <- lapply(reg_elements, function(i) 
		merge_records(ensembl_res[ensembl_res$regulatory_stable_id == i, ]))
	reg_vals_df <- Reduce(rbind, reg_vals)
	out <- data.frame(ensembl_reg_id = paste(reg_vals_df$ensembl_reg_id, collapse = "///"),
					  ensembl_reg_coordinates = paste(reg_vals_df$ensembl_reg_coordinates, collapse = "///"),
					  ensembl_reg_type = paste(reg_vals_df$ensembl_reg_type, collapse = "///"),
					  ensembl_reg_tissues = paste(reg_vals_df$ensembl_reg_tissues, collapse = "///"))
	out
}

#' Merge records for the same ENSEMBL regulatory element
#' 
#' This function collapses the activity status of a given an ENSEMBL regulatory 
#' element in different tissues. Notice that tissues identified as inactive will
#' not be reported.
#' 
#' @param tab Results from `biomaRt::getBM` for the same regulatory element
#' @return `data.frame` of one row after collapsing the 
merge_records <- function(tab){
	
	vec <- tab[1, , drop = FALSE]
	out <- data.frame(ensembl_reg_id = tab$regulatory_stable_id[1], 
					  ensembl_reg_coordinates = paste0(vec$chromosome_name, ":",
					  					 vec$chromosome_start, "-",
					  					 vec$chromosome_end), 
					  ensembl_reg_type = vec$feature_type_name)
	
	## Select activity and tissue
	subtab <- tab[, c("activity", "epigenome_name")]
	
	## Remove tissues without activity
	subtab <- subtab[!is.na(subtab$activity), ]
	states <- unique(subtab$activity)
	
	## Remove inactive from states
	states <- states[states != "INACTIVE"]
	
	state_vec <- lapply(states, function(x) {
		
		mini <- subtab[subtab$activity == x, ]
		paste(x, ":", paste(mini$epigenome_name, collapse = ";"))
	})
	
	out$ensembl_reg_tissues <- paste(unlist(state_vec), collapse = "/")
	out
}
