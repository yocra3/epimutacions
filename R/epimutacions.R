#' @export
epimutacions <- function(methy, method, min_cpg = 3, verbose = TRUE) #, num.cpgs = 10, pValue.cutoff = 0.01, 
					   #cutoff =0.1, outlier.score = 0.5, 
					   #nsamp = "deterministic",method = "manova")
{
	
	window_sz = 10
	offset_mean = 0.15
	offset_abs = 0.1
	
	bump_cutoff <-  0.1
	nsamp <- "exact"
	
	
	# Identify type of input and extract required data:
	#	* betas
	#	* sample's classification
	#	* feature annotation
	if(class(methy) == "GenomicRatioSet") {
		if(verbose) message("Input of type 'GeomicRatioSet")
		betas <- minfi::getBeta(methy)
		pd <- as.data.frame(SummarizedExperiment::colData(methy))
		fd <- as.data.frame(SummarizedExperiment::rowRanges(methy))
		rownames(fd) <- rownames(dta)
	} else if(class(methy) == "ExpressionSet") {
		if(verbose) message("Input of type 'ExpressionSet")
		betas <- Biobase::exprs(methy)
		pd <- Biobase::pData(methy)
		fd <- Biobase::fData(methy)
	} else {
		stop("Input data 'methy' must be a GenomicRatioSet or an ExpressionSet")
	}
	if(!"status" %in% colnames(pd)) {
		stop("epimutacions function required of column 'status' in sample's description with case/control labels")
	}
	
	# Identify the method to be used
	method <- charmatch(method, c("manova", "mlm", "iso.forest", "mahdist.mcd", "barbosa"))
	method <- c("manova", "mlm", "iso.forest", "mahdist.mcd", "barbosa")[method]
	if(method %in% c("iso.forest")) {
		stop("Method not implemented yet")
	}
	
	# Identify cases and controls
	ctr_sam <- rownames(pd)[pd$status == "control"]
	cas_sam <- rownames(pd)[pd$status == "case"]
	
	# Differentiate between methods that required region detection that the ones
	# that finds outliers to identify regions
	if(method %in% c("manova", "mlm", "mahdist.mcd")) {
		if(verbose) message(paste0("Selected method '", method, "' required of 'bumphunter'"))
		
		rst <- lapply(cas_sam, function(case) {
			# Prepare model to be evaluated
			pdi <- pd[ , "status", drop = FALSE]
			pdi$case <- rownames(pdi) == case
			model <- stats::model.matrix(~case, pdi)
			rm(pdi)
			
			# Run bumphunter for region partitioning
			bumps <- bumphunter::bumphunter(object = betas, design = model,
				pos = fd$start, chr = fd$seqnames, cutoff = bump_cutoff)$table
			bumps <- bumps[bumps$L >= min_cpg, ]
			bumps$sz <- bumps$end - bumps$start
			bumps <- bumps[bumps$sz < length(ctr_sam), ] # <--------------- TODO
			if(verbose) message(paste0(nrow(bumps), " candidate regions were found for case sample '", case, "'"))
			
			# Identify outliers according to selected method
			bump_out <- do.call(rbind, lapply(seq_len(nrow(bumps)), function(ii) {
				bump <- bumps[ii, ]
				beta_bump <- betas_from_bump(bump, fd, betas)
				
				if(method == "mahdist.mcd") {
					dst <- epi_mahdistmcd(beta_bump, nsamp)
					threshold <- sqrt(qchisq(p = 0.975, df = ncol(beta_bump)))
					
					# Find outliers according to the robust distance
					outliers <- which(dst$statistic >= threshold)
					
					# Report if there are ourliers
					bump$outlier <- case %in% outliers
					bump$outlier_score <- NA
					bump$outlier_significance <- NA
					bump$outlier_direction <- NA
					bump$CpG_ids <- paste(rownames(beta_bump), collapse = ",", sep = "")
					bump$sample <- case
					bump <- bump[bump$outlier, ]
					bump[ , c("chr", "start", "end", "sz", "L", "CpG_ids", "outlier_score", "outlier_significance", "outlier_direction", "sample")]
				} else if(method == "mlm") {
					sts <- epi_mlm(beta_bump, model)
					
					bump$outlier_score <- paste0(sts[1], "/", sts[2])
					bump$outlier_significance <- sts[3]
					bump$outlier_direction <- ifelse(bump$value < 0, "hypomethylation", "hypermethylation")
					bump$CpG_ids <- paste(rownames(beta_bump), collapse = ",", sep = "")
					bump$sample <- case
					bump[ , c("chr", "start", "end", "sz", "L", "CpG_ids", "outlier_score", "outlier_significance", "outlier_direction", "sample")]
				} else if(method == "manova") {
					sts <- epi_manova(beta_bump, model, case)
					
					bump$outlier_score <- paste0(sts[[1]][1], "/", sts[[1]][2])
					bump$outlier_significance <- sts[[1]][3]
					bump$outlier_direction <- ifelse(bump$value < 0, "hypomethylation", "hypermethylation")
					bump$CpG_ids <- paste(rownames(beta_bump), collapse = ",", sep = "")
					bump$sample <- case
					bump[ , c("chr", "start", "end", "sz", "L", "CpG_ids", "outlier_score", "outlier_significance", "outlier_direction", "sample")]
				}
			}))
		})
		
		return(rst)
	} else if(method == "barbosa") {
		# Compute reference statistics
		if(verbose) message("Calculating statistics from reference distribution required by Barbosa et. al. 2019")
		bctr_min <- apply(betas[ , ctr_sam], 1, min, na.rm = TRUE)
		bctr_max <- apply(betas[ , ctr_sam], 1, max, na.rm = TRUE)
		bctr_mean <- apply(betas[ , ctr_sam], 1, mean, na.rm = TRUE)
		bctr_prc <- apply(betas[ , ctr_sam], 1, quantile, probs = c(0.01, 0.99), na.rm = TRUE)
		bctr_pmin <- bctr_prc[1, ]
		bctr_pmax <- bctr_prc[2, ]
		rm(bctr_prc)
		#case <- betas[ , cas_sam[1], drop=FALSE]
		
		# Run region detection
		rst <- do.call(rbind, lapply(cas_sam, function(case) {
			x <- epi_barbosa(betas[ , case, drop=FALSE], fd, bctr_min, bctr_max, bctr_mean, 
				bctr_pmin, bctr_pmax, window_sz, min_cpg, offset_mean, offset_abs)
			x$sample <- case
			x
		}))
		
		return(rst)
	}
}