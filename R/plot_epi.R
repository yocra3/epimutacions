#' @export
plot_epi <- function(epi_res, method, sam, chr, genome = "hg19") {
	if(length(chr) != 1) {
		stop("Argument 'chr' must be of length one (aka. 'chr10')")
	}
	if(length(sam) != 1 | !sam %in% epi_res$sample) {
		stop("Argument 'sam' must be of length one and present in 'epi_res'")
	}
	epi_res <- epi_res[epi_res$sample == sam & epi_res$chromosome == chr, ]
	
	ideoTrack <- Gviz::IdeogramTrack(genome = genome, chromosome = chr)
	gatrack <- Gviz::GenomeAxisTrack()
	atrack <- Gviz::AnnotationTrack(epi_res, name = "Epimutations")
	
	if(method %in% c("barbosa", "mahdistmcd")) {
		Gviz::plotTracks(list(ideoTrack, gatrack, atrack, qtrack))	
	} else if(method %in% c("manova", "mlm")) {
		qtrack = Gviz::DataTrack(epi_res, data = -log10(epi_res$outlier_significance),
								 name = "Significance [-log10(pval)]", genome = genome)
		Gviz::plotTracks(list(ideoTrack, gatrack, atrack, qtrack))
	} else if(method == "isoforest") {
		qtrack = Gviz::DataTrack(epi_res, data = epi_res$outlier_score,
								 name = "Score", genome = genome)
		Gviz::plotTracks(list(ideoTrack, gatrack, atrack, qtrack))
	}
}