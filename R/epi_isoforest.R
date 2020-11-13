epi_isoforest <- function(mixture, case_id) {
	mixture <- t(mixture)
	#Generate train and test(sample with suspected disease) data frame
	train <- mixture[row.names(mixture) != case_id, ]
	test <- mixture[case_id, , drop = FALSE]
	
	#Run the isolation forest methods 
	iso <- isotree::isolation.forest(train)
	
	#Predict
	score <- predict(iso, test)
	
	return(score)
}

res_isoforest <- function(bump, beta_bump, sts, case) {
	bump$outlier_score <- sts
	bump$outlier_significance <- NA
	bump$outlier_direction <- ifelse(bump$value < 0, "hypomethylation", "hypermethylation")
	bump$CpG_ids <- paste(rownames(beta_bump), collapse = ",", sep = "")
	bump$sample <- case
	bump[ , c("chr", "start", "end", "sz", "L", "CpG_ids", "outlier_score", "outlier_significance", "outlier_direction", "sample")]
}