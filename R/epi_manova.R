epi_manova <-  function(mixture, model, case_id){
	mixture <- t(mixture)
	mod <- manova(mixture ~ model)
	mod_summary <- summary(mod)$stats
	statistics <- mod_summary[1, c("approx F", "Pillai","Pr(>F)")]
	
	# Calculate the beta mean difference
	keep_case <- which(case_id %in% rownames(mixture))
	case <- mixture[keep_case,]
	controls <- mixture[-keep_case, ]
	coltrols_mean <- colMeans(controls)
	beta_mean_difference <- coltrols_mean - case
	
	output <- list(statistics, beta_mean_difference)
	return(output)
}

res_manova <- function(bump, beta_bump, sts, case) {
	bump$outlier_score <- paste0(sts[[1]][1], "/", sts[[1]][2])
	bump$outlier_significance <- sts[[1]][3]
	bump$outlier_direction <- ifelse(bump$value < 0, "hypomethylation", "hypermethylation")
	bump$CpG_ids <- paste(rownames(beta_bump), collapse = ",", sep = "")
	bump$sample <- case
	bump[ , c("chr", "start", "end", "sz", "L", "CpG_ids", "outlier_score", "outlier_significance", "outlier_direction", "sample")]
}