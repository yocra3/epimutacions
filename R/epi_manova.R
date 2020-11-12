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