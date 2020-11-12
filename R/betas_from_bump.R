betas_from_bump <- function(bump, fd, betas) {
	cpgs <- rownames(fd[fd$seqnames == bump$chr & fd$start >= bump$start & fd$start <= bump$end, ])
	return(betas[cpgs, ])
}
