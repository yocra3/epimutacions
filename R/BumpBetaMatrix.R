#' Function to extract the methylation of a specific region
#' 
#' This function extracts the methylation from an GenomicRatioSet or
#' ExpressionSet object given a region defined as a bumps object.
#' 
#' @param bumps The region to be extracted as a bumps object.
#' @param dataset The full methylation set as GenomicRatioSet or ExpressionSet.
#' @return The betas as a matrix.
#' @export
get_betas <- function(bump, dataset) {
	# Identify the class of the metylation data set
	data_type <- charmatch(class(dataset), c("GenomicRatioSet", "ExpressionSet"))
	if(is.na(data_type)) {
		stop("Argument 'dataset' must be a 'GenomicRatioSet' or an 'ExpressionSet'")
	}
	# Extract the beta values
	if(data_type == 1) {
		betas <- minfi::getBeta(dataset[bump$indexStart:bump$indexEnd, ]) 
	} else {
		betas <- Biobase::exprs(dataset[bump$indexStart:bump$indexEnd, ])
	}
	# Remove NAs
	betas <- na.omit(betas)
	
	# Return transposed matrix
	return(t(betas))
}

#' Return the average regional methylation difference between a case and controls.
#'
#' @param betas (matrix) The beta values. Columns are CpGs, rows are samples.
#' @param case_id (string) The rowname in betas for the case of interest.
#'
#' @return (numeric) mean_beta(case) - mean_beta(controls)
#' @keywords internal
ave_beta_case_minus_controls <- function(betas, case_id){
	case_row <- which(rownames(betas) %in% case_id)
	average_betas_per_sample <- rowMeans(betas)
	return(
		average_betas_per_sample[case_row] - 
			mean(average_betas_per_sample[-case_row])
	)
}