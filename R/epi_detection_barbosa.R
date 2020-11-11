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
	if(ncol(cases) != 1) {
		stop("Epimutation detection with berbosa approach can only works with a singe proband")
	}
	
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
		flag_q_sup = bcas[, 1] >= bctr_pmax & bcas[, 1] >= bctr_max + offset_abs,
		flag_q_inf = bcas[, 1] <= bctr_pmin & bcas[, 1] <= bctr_min - offset_abs,
		flag_m_sup = bcas[, 1] >= bctr_mean + offset_mean,
		flag_m_inf = bcas[, 1] <= bctr_mean - offset_mean,
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
		# Order the input first by chromosome and then by position
		flag_df <- flag_df[with(flag_df, order(chr, pos)), ]
		
		x <- do.call(rbind, lapply(unique(flag_df$chr), function(chr) {
			# Subset by chromosome
			red_df <- flag_df[flag_df$chr == chr, ]
			if(nrow(red_df) < N) {
				return(data.frame(chr=NA, pos=NA, probands=NA, region=NA))
			}
			
			# Get the position of the previous and next probe for each probe in the
			# data.frame. The first and last position get its own position minus/plus
			# the window size to be sure to include them in the resulting data.frame.
			red_df$pos_next <- c(red_df$pos[seq(2, nrow(red_df))], red_df$pos[nrow(red_df)] + window_sz + 1)
			red_df$pos_prev <- c(red_df$pos[1] - window_sz - 1, red_df$pos[seq(1, nrow(red_df) - 1)])
			
			# We add two columns indicating if a probe is within the window size
			# range with its previous and with its next probe
			red_df$in_prev <- red_df$pos - red_df$pos_prev <= window_sz
			red_df$in_next <- red_df$pos_next - red_df$pos <= window_sz
			
			# We drop all the probes that do not have a previous nor a next
			# probe within the range of interest
			red_df <- red_df[red_df$in_prev | red_df$in_next, ]
			if(nrow(red_df) < N) {
				return(data.frame(chr=NA, pos=NA, probands=NA, region=NA))
			}
			
			# Using the cumsum function and by negating the content of the "in_next"
			# column we can define the regions of CpGs within the range since they
			# will be tagged with the same number
			red_df$cum <- cumsum(!red_df$in_next)
			
			# correct the base position of the change in the region
			red_df$cum2 <- red_df$cum
			for(ii in seq(2, nrow(red_df))) {
				if(red_df$cum[ii] != red_df$cum[ii - 1] & red_df$in_prev[ii] & !red_df$in_next[ii]) {
					red_df$cum2[ii] <- red_df$cum[ii] - 1
				}
			}
			
			# Computing the frequency of each "number" assign to the region we can 
			# know how may probes are in it. We can use this frequency to filter out
			# those regions with less probes than given by N.
			# We also give to the regions a proper name.
			fr <- data.frame(table(red_df$cum2), stringsAsFactors = FALSE)
			fr <- as.numeric(as.character(fr$Var1[fr$Freq >= N]))
			if(length(fr) > 0) {
				fr <- data.frame(current = fr, new = paste0(pref, "_", chr, "_", seq_len(length(fr))))
				red_df <- red_df[red_df$cum2 %in% fr$current, ]
				rownames(fr) <- paste0("O", fr$current)
				red_df$region <- fr[paste0("O", red_df$cum2), "new"]
				
				
				# Since the first and last probe in a chromosome will have TRUE in
				# prev or next distance we need to be sure to drop them if they
				# are not in the window
				red_df$dist_next <- red_df$pos_next - red_df$pos
				red_df$dist_next[length(red_df$dist_next)] <- 0
				red_df <- red_df[red_df$dist_next <= window_sz, ]
				
				# We drop the columns with the flags used for the outlier and region
				# detection
				red_df <- red_df[ , c("chr", "pos", "probands", "region")]
				return(red_df)
			} else {
				return(data.frame(chr=NA, pos=NA, probands=NA, region=NA))
			}
		}))
		
		return(x[!is.na(x$chr), ])
	}
	
	# We select the CpGs according to the direction of the outlier
	flag_sup <- flag_result[flag_result$flag_qm_sup, ]
	flag_inf <- flag_result[flag_result$flag_qm_inf, ]
	
	# We identify the regions taking into account the direction
	reg_sup <- get_regions(flag_sup, pref = "Rs")
	reg_inf <- get_regions(flag_inf, pref = "Ri")
	
	# We add a column indicating the direction of the regions/outliers
	reg_sup$outlier_score <- "hypermethylation"
	reg_sup$CpG_ids <- rownames(reg_sup)
	reg_inf$outlier_score <- "hypomethylation"
	reg_inf$CpG_ids <- rownames(reg_inf)
	
	collapse_regions <- function(flag_df) {
		do.call(rbind, lapply(unique(flag_df$region), function(reg) {
			x <- flag_df[flag_df$region == reg, ]
			data.frame(
				chr = x$chr[1],
				start = min(x$pos),
				end = max(x$pos),
				length = max(x$pos) - min(x$pos),
				N_CpGs = nrow(x),
				CpG_ids = paste(x$CpG_ids, collapse = "/", sep = ""),
				outlier_method = "barbosa",
				outlier_score = x$outlier_score[1]
			)
				
		}))
	}
	
	# We collapse the CpGs in regions and format the output
	clean_sup <- collapse_regions(reg_sup)
	clean_inf <- collapse_regions(reg_inf)
	
	return(rbind(clean_inf, clean_sup))
}



