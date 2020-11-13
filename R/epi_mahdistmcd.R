#' Epimutation detection based on Outlier Detection with Robust Mahalanobis
#' 
#' The implementation of the method here is based on the discusion in this
#' thred of StakOverflow:
#' https://stats.stackexchange.com/questions/137263/multivariate-outlier-detection-with-robust-mahalanobis
epi_mahdistmcd <- function(mixture, nsamp = c("best", "exact", "deterministic")) {
	nsamp <- charmatch(nsamp, c("best", "exact", "deterministic"))
	nsamp <- c("best", "exact", "deterministic")[nsamp]
	if(is.na(nsamp)) {
		stop("Argument 'nsamp' shuld be 'best', 'exact', 'deterministic'")
	}
	
	# Transpose input methylation beta matrix to have the samples as rows and 
	# the CpGs as columns
	mixture <- t(mixture)
	
	# Run the MCD model on the transposed matrix
	mcd <- robustbase::covMcd(mixture, nsamp = nsamp)
	
	# Get MCD stimate of location
	mean_mcd <- mcd$center
	
	# Get MCD estimate scatter
	cov_mcd <- mcd$cov
	
	# Get inverse of scatter
	cov_mcd_inv <- solve(cov_mcd)
	
	# Compute the robust distance between samples
	robust_dist <- apply(mixture, 1, function(x){
		x <- (x - mean_mcd)
		dist <- sqrt((t(x)  %*% cov_mcd_inv %*% x))
		return(dist)
	})
	
	return(data.frame(ID = names(robust_dist), statistic = robust_dist))
}


res_mahdistmcd <- function(case, bump, beta_bump, outliers) {
	bump$outlier <- case %in% outliers
	bump$outlier_score <- NA
	bump$outlier_significance <- NA
	bump$outlier_direction <- NA
	bump$CpG_ids <- paste(rownames(beta_bump), collapse = ",", sep = "")
	bump$sample <- case
	bump <- bump[bump$outlier, ]
	bump[ , c("chr", "start", "end", "sz", "L", "CpG_ids", "outlier_score", "outlier_significance", "outlier_direction", "sample")]
}