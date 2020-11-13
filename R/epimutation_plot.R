#' Function to plot selected CpGs from epimutacions package. NOT WORKING YET!
#' 
#' This function plots a methylated region in genome coordinates using reshape package
#' @param data getRationSetObject.
#' @param res results from epimutacions package.
#' @param samp_ID column name of sample ids . Default: 'samp_ID'.
#' @param cpg_col column name of sample ids . Default: 'CpG_ids'.
#' @return List of plots.
#' @export
#' 
#' # 2 inputs, 1 the original and the other one from the methods with the signficant ones, then selection of what CpG we want to see
# add chromosome
function(data, res, samp_ID='samp_ID', cpg_col="CpG_ids"){
	# FUNCTION NOT FULLY TESTED. TO TRY IT OUT
	# MANUAL LOAD USING: data(genomicratioset) and then res <- epimutations(genomicratioset)
	# get beta values and its correspondent genomic coordinates
	beta_values <-  minfi::getBeta(genomicratioset) 
	beta_coord <- SummarizedExperiment::rowRanges(genomicratioset)
	
	# separate controls from case
	beta_values_control <- beta_values[-which(colnames(beta_values)%in%res[[samp_ID]]),]
	cases <- which(colnames(beta_values)%in%res[[samp_ID]])
	
	# if cases is only 1 need to be processed differently to make sure we have a data frame
	beta_values_cases <- beta_values[,cases, drop=F]
	
	# make a plot for each case and store it in a list
	plots <- list()
	for (case in colnames(beta_values)[cases]){
		# make plot for each coordinate
		subplots <- list()
		
		processed_case <- res[res[[samp_ID]] == case]
		 
		selected_islands <- beta_values[processed,]
		apply(X = processed_case, MARGIN = 1, FUN = function(row){
			
			processed_cpgs <- unlist(strsplit(row[[cpg_col]], ','))
			selected_islands <- beta_values[processed_cpgs,case, drop=F]
			selected_coord <- beta_coord[processed_cpgs,]
			
			
			data_controls <- cbind(data.frame(coord=data.frame(ranges(selected_coord))$start), beta_values_control[processed_cpgs,])
			data_case <- cbind(data.frame(coord=data.frame(ranges(selected_coord))$start), selected_islands[,case, drop=F])
			colnames(data_case) <- c('coord', 'case')
			plot(runif(10), runif(10),
				 xlim = c(row$start, row$end), 
				 ylim=c(0,1),
				 type="n",  #hide the points
				 ylab="Methylation Fraction", xlab="Genomic coordinates")
			apply(X = data_controls, MARGIN = 2, FUN = function(x){
				lines(x=data_controls$coord, y=x, type='o', col="#9d95a1CC")
			})
			lines(x=data_case$coord, y=data_case$case, type='o', pch=19, col='red')
			
			
			
		})
	}
	
	}

