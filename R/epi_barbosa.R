#' Method implementing Barbosa et. al. 2019 to detect epimutations
#' 
#' The implementation of the method here is based on teh description of the
#' method used in:
#' Barbosa M, Joshi RS, Garg P, et al. Identification of rare de novo epigenetic
#'  variations in congenital disorders. Nat Commun. 2018;9(1):2064. Published 
#'  2018 May 25. doi:10.1038/s41467-018-04540-x
epi_barbosa <- function(case, fd, bctr_min, bctr_max, bctr_mean, bctr_pmin, 
	bctr_pmax, window_sz = 1000, N = 3, offset_mean = 0.15, offset_abs = 0.1) {
	# Check that there is a single proband
	if(ncol(case) != 1) {
		stop("Epimutation detection with barbosa approach can only works with a singe proband")
	}
	
	# Identify outlier at the proband side using reference statistics
	flag_result <- data.frame(
		chr = as.character(fd[rownames(case), "seqnames"]),
		pos = fd[rownames(case), "start"],
		flag_q_sup = case[, 1] >= bctr_pmax & case[, 1] >= bctr_max + offset_abs,
		flag_q_inf = case[, 1] <= bctr_pmin & case[, 1] <= bctr_min - offset_abs,
		flag_m_sup = case[, 1] >= bctr_mean + offset_mean,
		flag_m_inf = case[, 1] <= bctr_mean - offset_mean,
		stringsAsFactors = FALSE
	)
	flag_result$flag_q_sup[is.na(flag_result$flag_q_sup)] <- FALSE
	flag_result$flag_m_sup[is.na(flag_result$flag_m_sup)] <- FALSE
	flag_result$flag_q_inf[is.na(flag_result$flag_q_inf)] <- FALSE
	flag_result$flag_m_inf[is.na(flag_result$flag_m_inf)] <- FALSE
	flag_result$flag_qm_sup <- flag_result$flag_q_sup & flag_result$flag_m_sup
	flag_result$flag_qm_inf <- flag_result$flag_q_inf & flag_result$flag_m_inf
	
	# Function used to detect regions of N CpGs closer than window size
	get_regions <- function(flag_df, window_sz = 1000, N = 3, pref = "R") {
		if(nrow(flag_df) < N) {
			return(data.frame(chr=NA, pos=NA, region=NA))
		}
		
		# Order the input first by chromosome and then by position
		flag_df <- flag_df[with(flag_df, order(chr, pos)), ]
		
		x <- do.call(rbind, lapply(unique(flag_df$chr), function(chr) {
			# Subset by chromosome
			red_df <- flag_df[flag_df$chr == chr, ]
			if(nrow(red_df) < N) {
				return(data.frame(chr=NA, pos=NA, region=NA))
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
				return(data.frame(chr=NA, pos=NA, region=NA))
			}
			
			# Using the cumsum function and by negating the content of the "in_next"
			# column we can define the regions of CpGs within the range since they
			# will be tagged with the same number
			red_df$cum <- cumsum(!red_df$in_next)
			
			# Correct the base position of the change in the region
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
				red_df <- red_df[ , c("chr", "pos", "region")]
				return(red_df)
			} else {
				return(data.frame(chr=NA, pos=NA, region=NA))
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
	reg_sup$outlier_direction <- "hypermethylation"
	reg_sup$CpG_ids <- rownames(reg_sup)
	reg_inf$outlier_direction <- "hypomethylation"
	reg_inf$CpG_ids <- rownames(reg_inf)
	
	collapse_regions <- function(flag_df) {
		do.call(rbind, lapply(unique(flag_df$region), function(reg) {
			x <- flag_df[flag_df$region == reg, ]
			data.frame(
				chromosome = x$chr[1],
				start = min(x$pos),
				end = max(x$pos),
				length = max(x$pos) - min(x$pos),
				N_CpGs = nrow(x),
				CpG_ids = paste(x$CpG_ids, collapse = ",", sep = ""),
				outlier_score = NA,
				outlier_significance = NA,
				outlier_direction = x$outlier_direction[1]
			)
			
		}))
	}
	
	# We collapse the CpGs in regions and format the output
	clean_sup <- collapse_regions(reg_sup)
	clean_sup <- clean_sup[!is.na(clean_sup$chromosome), ]
	clean_inf <- collapse_regions(reg_inf)
	clean_inf <- clean_inf[!is.na(clean_inf$chromosome), ]
	
	return(rbind(clean_inf, clean_sup))
}



