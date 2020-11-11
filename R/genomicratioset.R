#' \code{GenomicRatioSet} class object, created from
#' methylation_matrix dataset
#'
#' \code{GenomicRatioSet} object, annotated with
#' \code{IlluminaHumanMethylation450kanno.ilmn12.hg19} data
#' with contiguous CpG choosing on 15 chr. 
#' 
#' 
#' Dimensions: The dataset consist of 20 rows (representing CpG islands)
#' and 5 columns (representing samples)
#' 
#' Sample choosing was done with subseting the default Locations
#' dataframe:
#' \code{chr15 <- data.frame(Locations) %>%
#' filter(chr == 'chr15')}
#' After that, new dataframe was ordered by position in ascending order
#' and slice with 20 CpG names was taken:
#' \code{row_n <- rownames(arrange(chr15, chr15$pos)[10000:10019,])}
#' 
#' The geneticratioset was created via:
#' \code{geneticratioset <- makeGenomicRatioSetFromMatrix(methylation_matrix, 
#' rownames = row_n, pData = data.frame(id = c("case_1", "case_2", "case_3", 
#' "case_4", "case_5")))}
#'
#' @usage data("genomicratioset")
#' @return A \code{GenomicRatioSet} object.
#' @examples
#' data("methylation_matrix")
#' dim(methylation_matrix)
#' sampleNames(genomicratioset)
#'
"genomicratioset"