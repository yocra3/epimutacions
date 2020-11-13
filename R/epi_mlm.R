#' Epimutation detection using Multivariate linear models
epi_mlm <- function(mixture, model, case) {
	mod <- mlm::mlm(t(mixture) ~ model[,2])
	statistics <- mod$aov.tab[1, c("F value", "R2", "Pr(>F)")]
	return(statistics)
}

res_mlm <- function(bump, beta_bump, sts, case) {
	bump$outlier_score <- paste0(sts[1], "/", sts[2])
	bump$outlier_significance <- sts[3]
	bump$outlier_direction <- ifelse(bump$value < 0, "hypomethylation", "hypermethylation")
	bump$CpG_ids <- paste(rownames(beta_bump), collapse = ",", sep = "")
	bump$sample <- case
	bump[ , c("chr", "start", "end", "sz", "L", "CpG_ids", "outlier_score", "outlier_significance", "outlier_direction", "sample")]
}