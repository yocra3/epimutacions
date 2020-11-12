#' Epimutation detection using Multivariate linear models
#' 
epi_mlm <- function(mixture, model) {
	mod <- mlm::mlm(t(mixture) ~ model[,2])
	statistics <- mod$aov.tab[1, c("F value", "R2", "Pr(>F)")]
	return(statistics)
}