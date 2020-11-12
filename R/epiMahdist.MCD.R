#' Epimutation detection based on Outlier Detection with Robust Mahalanobis
#' 
#' The implementation of the method here is based on the discusion in this
#' thred of StakOverflow:
#' https://stats.stackexchange.com/questions/137263/multivariate-outlier-detection-with-robust-mahalanobis
#' 
#' @param betas matrix A matrix with the beta values of controls with a single 
#' case (proband) to be tested for outlier with CpGs as rows and samples 
#' (controls + 1 case ) as columns.
#' @param nsamp character/integer Number of subsets used for initial estimates.
#' It can take a stratgey ("best", "exact", or "deterministic") or an exact
#' number.
#' @export
epiMahdist.MCD <- function(betas, nsamp = c("best", "exact", "deterministic")) {
	if(is.null(betas)) {
		stop("'betas' data frame must be introduced")
	}
	
	val <- match(nsamp, c("best", "exact", "deterministic"))
	if(length(val) != 1) {
		stop("nsamp shuld be 'best', 'exact', 'deterministic'")
	}
	nsamp <- c("best", "exact", "deterministic")[val]
	
	if(is.null(sample) | length(sample)!= 1) {
		stop("One 'sample' with suspected rare disease must be introduced")
	}
	
	# Transpose input methylation beta matrix to have the samples as rows and 
	# the CpGs as columns
	betas_t <- t(betas)
	
	# Run the MCD model on the transposed matrix
	mcd <- robustbase::covMcd(betas_t, nsamp = nsamp)
	
	# Get MCD stimate of location
	mean_mcd <- mcd$center
	
	# Get MCD estimate scatter
	cov_mcd <- mcd$cov
	
	# Get inverse of scatter
	cov_mcd_inv <- solve(cov_mcd)
	
	# Compute the robust distance between samples
	robust_dist <- apply(betas_t, 1, function(x){
		x <- (x - mean_mcd)
		dist <- sqrt((t(x)  %*% cov_mcd_inv %*% x))
		return(dist)
	})
	
	return(data.frame(ID = names(robust_dist), statistic = robust_dist))
}




# #set cutoff using chi square distribution
# threshold <- sqrt(qchisq(p = 0.975, df = ncol(betas)))
# 
# #find outliers
# outliers <-  which(robust_dist >= threshold)
# 
# #suspected rare disease sample within the outliers
# 
# rd.sample.out <- sample %in% names(outliers)
# 
# return(rd.sample.out)