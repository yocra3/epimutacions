#' MEthod implementing Barbosa et. al. 2019 to detect epimutations
#' 
#' @param cases GenomicRatioSet with cases (probands) methylation data
#' @param controls GenomicRatioSet with controls methylation data to be used as
#' reference.
#' @param window integer defining the window size where N CpGs have to be to 
#' define a region (default: 1000)
#' @param N integerdefining the number of CpGs to closer than window bases to
#' define a region (default: 3)
#' @export
epi_detection_barbosa <- function(cases, controls, window = 1000) {
	# Compute requires statistics from controls for each probe:
	#    * min reference beta value
	#    * max reference beta value
	#    * 99th percentile of reference beta value
	#    * 0.1st percentile of reference beta value
	#    * mean of reference beta value
	bctr <- minfi::getBeta(controls)
	bctr_min <- apply(bctr, 1, min, na.rm = TRUE)
	bctr_max <- apply(bctr, 1, max, na.rm = TRUE)
	bctr_mean <- apply(bctr, 1, mean, na.rm = TRUE)
	bctr_pmin <- apply(bctr, 1, quantile, probs = 0.01, na.rm = TRUE)
	bctr_pmax <- apply(bctr, 1, quantile, probs = 0.99, na.rm = TRUE)
	rm(bctr)
	
	# Identify outlier at the probands side using previous statistics
	bcas <- minfi::getBeta(cases)
	flag_result <- data.frame(
		chr = as.character(seqnames(rowRanges(cases))),
		pos = start(rowRanges(cases)),
		probands = bcas[ , 1],
		flag_q_sup = bcas[, 1] >= bctr_pmax & bcas[, 1] >= bctr_max,
		flag_q_inf = bcas[, 1] <= bctr_pmin & bcas[, 1] <= bctr_min,
		flag_m_sup = bcas[, 1] >= bctr_mean,
		flag_m_inf = bcas[, 1] <= bctr_mean,
		stringsAsFactors = FALSE
	)
	flag_result$flag_q_sup[is.na(flag_result$flag_q_sup)] <- FALSE
	flag_result$flag_m_sup[is.na(flag_result$flag_m_sup)] <- FALSE
	flag_result$flag_q_inf[is.na(flag_result$flag_q_inf)] <- FALSE
	flag_result$flag_m_inf[is.na(flag_result$flag_m_inf)] <- FALSE
	flag_result$flag_qm_sup <- flag_result$flag_q_sup & flag_result$flag_m_sup
	flag_result$flag_qm_inf <- flag_result$flag_q_inf & flag_result$flag_m_inf
	rm(bctr_min, bctr_max, bctr_mean, bctr_pmin, bctr_pmax)
	rm(bcas)
	
	# Function used to detect regions of N CpGs closer than window size
	get_regions <- function(flag_df, N = 3) {
		regions <- flag_df[1, ]
		regions$region <- 1
		r_cnt <- 1
		r_ii <- 1
		for(ii in seq(2, nrow(flag_df))) {
			# if the new CpG is from another chromosome, we create a new region (r_cnt),
			# we add the new CpG, and we increse the row counter (r_ii)
			if(regions$chr[r_ii] != flag_df$chr[ii]) {
				r_cnt <- r_cnt + 1
				regions <- rbind(regions, data.frame(flag_df[ii, ], region = r_cnt))
				r_ii <- r_ii + 1
			} else {
				# if the new CpG is in the window taking into account the last include
				# GpG, we include it and increase the row-counter (r_ii)
				if(regions$pos[r_ii] + window >= flag_df$pos[ii]) {
					regions <- rbind(regions, data.frame(flag_df[ii, ], region = r_cnt))
					r_ii <- r_ii + 1
				}
				# if the new CpG is out the window taking into account the last include
				# GpG, we include the region counter (r_cnt) it and increase the 
				# row-counter (r_ii)
				if(regions$pos[r_ii] + window < flag_df$pos[ii]) {
					r_cnt <- r_cnt + 1
					regions <- rbind(regions, data.frame(flag_df[ii, ], region = r_cnt))
					r_ii <- r_ii + 1
				}
			}
		}
		
		fr <- data.frame(table(regions$region), stringsAsFactors = FALSE)
		fr <- as.numeric(fr$Var1[fr$Freq > N])
		
		regions <- regions[regions$region %in% fr, ]
		regions$region <- paste0("R", regions$region)
		return(regions)
	}
	
	flag_sup <- flag_result[flag_result$flag_qm_sup, ]
	flag_inf <- flag_result[flag_result$flag_qm_inf, ]
	
	reg_sup <- get_regions(flag_sup)
	reg_inf <- get_regions(flag_inf)
	
	reg_sup$direction <- "+"
	reg_inf$direction <- "-"
	
	return(rbind(reg_inf, reg_sup))
}



