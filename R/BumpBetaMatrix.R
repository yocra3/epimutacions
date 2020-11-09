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
	if(type == 1) {
		betas <- minfi::getBeta(set[bump$indexStart:bump$indexEnd, ]) 
	} else {
		betas <- Biobase::exprs(set[bump$indexStart:bump$indexEnd, ])
	}
	# Remove NAs
	betas <- na.omit(betas)
	
	# Return transposed matrix
	return(t(betas))
}
