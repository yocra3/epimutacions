#' Method implementing Barbosa et. al. 2019 to detect epimutations
#' 
#' This function implements the method for epimutation detection by Berbose et.
#' al. 2019. The process consist in identifying regions of CpGs that are
#' hypermethylated or hypometilated as:
#' 
#'  1. Hypermethylation: the case presents, in a 1kb window, probes that 
#'     fulfill both of the following criteria:
#'       A. At least 3 probes that each have betas above the 99.9th percentile 
#'          of the reference distribution for that probe, and are >= 0.15 above 
#'          the reference mean. 
#'       B. At least 1 probe with a beta value >= 0.1 above the maximum observed
#'          in reference for that probe.
#'  2. Hypomethylation: the case presents, in a 1kb window, probes that 
#'     fulfill both of the following criteria:
#'       A. At least 3 probes that each have beta values below the 0.1th 
#'          percentile of the reference distribution for that probe, and are 
#'          >= 0.15 below the reference mean. 
#'       B. At least 1 probe with a β value >=0.1 below the minimum observed in
#'          reference for that probe.
#' 
#' @param cases GenomicRatioSet with cases (probands) methylation data
#' @param controls GenomicRatioSet with controls methylation data to be used as
#' reference.
#' @param window_sz integer defining the window size where N CpGs have to be to 
#' define a region (default: 1000)
#' @param N integer defining the number of CpGs to closer than window bases to
#' define a region (default: 3)
#' @param offset_mean float defines the offset for the mean comparison, 
#' understood as that a probe has to be offset_mean over the reference mean to 
#' be considered outlier
#' @param offset_abs float defining the offset for the min/max comparison,
#' understood as that a probe has to be offset_abs under/over the min/max
#' value of a probe in the reference to be considered outlier
#' @export
epi_detection_barbosa <- function(cases, controls, window_sz = 1000, N = 3, offset_mean = 0.15, offset_abs = 0.1) {
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
	get_regions <- function(flag_df, window_sz = 1000, N = 3, pref = "R") {
		
		flag_df$pos_next <- c(flag_df$pos[seq(2, nrow(flag_df))], flag_df$pos[nrow(flag_df)] + window_sz)
		flag_df$pos_prev <- c(flag_df$pos[1] - window_sz, flag_df$pos[seq(1, nrow(flag_df) - 1)])
		
		flag_df$in_prev <- flag_df$pos - flag_df$pos_prev <= window_sz
		flag_df$in_next <- flag_df$pos_next - flag_df$pos <= window_sz
		
		flag_df <- flag_df[flag_df$in_prev | flag_df$in_next, ]
		
		flag_df$cum <- cumsum(!flag_df$in_next)
		
		fr <- data.frame(table(flag_df$cum), stringsAsFactors = FALSE)
		fr <- as.numeric(fr$Var1[fr$Freq > N])
		fr <- data.frame(current = fr, new = paste0(pref, seq_len(length(fr))))
		flag_df <- flag_df[flag_df$cum %in% fr$current, ]
		rownames(fr) <- paste("O", fr$current)
		flag_df$region <- fr[paste("O", flag_df$cum), "new"]
		flag_df <- flag_df[ , c("chr", "pos", "probands", "region")]
		
		return(flag_df)
		
		# regions <- flag_df[1, ]
		# regions$region <- 1
		# r_cnt <- 1
		# r_ii <- 1
		# for(ii in seq(2, nrow(flag_df))) {
		# 	# if the new CpG is from another chromosome, we create a new region (r_cnt),
		# 	# we add the new CpG, and we increse the row counter (r_ii)
		# 	if(regions$chr[r_ii] != flag_df$chr[ii]) {
		# 		r_cnt <- r_cnt + 1
		# 		regions <- rbind(regions, data.frame(flag_df[ii, ], region = r_cnt))
		# 		r_ii <- r_ii + 1
		# 	} else {
		# 		# if the new CpG is in the window taking into account the last include
		# 		# GpG, we include it and increase the row-counter (r_ii)
		# 		if(regions$pos[r_ii] + window >= flag_df$pos[ii]) {
		# 			regions <- rbind(regions, data.frame(flag_df[ii, ], region = r_cnt))
		# 			r_ii <- r_ii + 1
		# 		}
		# 		# if the new CpG is out the window taking into account the last include
		# 		# GpG, we include the region counter (r_cnt), and increase the 
		# 		# row-counter (r_ii)
		# 		if(regions$pos[r_ii] + window < flag_df$pos[ii]) {
		# 			r_cnt <- r_cnt + 1
		# 			regions <- rbind(regions, data.frame(flag_df[ii, ], region = r_cnt))
		# 			r_ii <- r_ii + 1
		# 		}
		# 	}
		# }
		# 
		# # We remove all the regions that does not have N CpGs in them
		# fr <- data.frame(table(regions$region), stringsAsFactors = FALSE)
		# fr <- as.numeric(fr$Var1[fr$Freq > N])
		# regions <- regions[regions$region %in% fr, ]
		# 
		# # We rename the regions by adding "R"
		# regions$region <- paste0("R", regions$region)
		# return(regions)
	}
	
	# We select the CpGs according to the direction of the outlier
	flag_sup <- flag_result[flag_result$flag_qm_sup, ]
	flag_inf <- flag_result[flag_result$flag_qm_inf, ]
	
	# We identify the regions taking into account the direction
	reg_sup <- get_regions(flag_sup)
	reg_inf <- get_regions(flag_inf)
	
	# We add a column indicating the direction of the regions/outliers
	reg_sup$direction <- "+"
	reg_inf$direction <- "-"
	
	return(rbind(reg_inf, reg_sup))
}



