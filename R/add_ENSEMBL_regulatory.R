#' Add ENSEMBL regulatory regions to epimutations
#' 
#' @param epimutations Result from `EpiMutations`
#' @param build Build used to define epimutations coordinates. By default, it is 37,
#' corresponding to Illumina annotation.
add_ENSEMBL_regulatory <- function(epimutations, build = "37"){
	
	## Create connection to ENSEMBL 
	mart <- biomaRt::useEnsembl(biomart = "regulation", GRCh = build)
	ensembl <- biomaRt::useDataset(dataset = "hsapiens_regulatory_feature", 
								   mart = mart)

	reg_res <- lapply(seq_len(nrow(epimutations)), function(i) {
		get_ENSEMBL_data(epimutations[i, "chromosome"], 
						 epimutations[i, "start"],
						 epimutations[i, "end"], mart = mart)
	})
	reg_res_df <- Reduce(rbind, reg_res)
	cbind(epimutations, reg_res_df)

}

get_ENSEMBL_data <- function(chromosome, start, end, mart){
	bm <- biomaRt::getBM(attributes = c("activity", "regulatory_stable_id", 
								  "chromosome_name", "chromosome_start",
								  "chromosome_end", "feature_type_name",
								  "epigenome_name"), 
				   filters = c('chromosome_name','start','end'),
				   values = list(chromosome, start, end),
				   mart = ensembl)
	out_ens <- process_ENSEMBL_results(bm)
	out_ens
}

#' Process data from ENSEMBL to combine results from the same regulatory elements in 
#' a unique record
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
