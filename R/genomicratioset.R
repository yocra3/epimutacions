#' \code{GenomicRatioSet} class object, created from
#' methylation_matrix dataset
#'
#' \code{GenomicRatioSet} object, annotated with
#' \code{IlluminaHumanMethylation450kanno.ilmn12.hg19} data
#' with random sample choosing (via \code{sample(485512, 20)}). 
#' Only rows for CpG_4, _5, _6, _7, _8 are chosen
#' with \code{rownames(Locations)[3000:3004]} in order to simulate
#' real hit with continous CpG island methylation difference.
#' 
#' Dimensions: The dataset consist of 20 rows (representing CpG islands)
#' and 5 columns (representing samples)
#'
#' @usage data("genomicratioset")
#' @return A \code{GenomicRatioSet} object.
#' @examples
#' data("methylation_matrix")
#' dim(methylation_matrix)
#' sampleNames(genomicratioset)
#'
"genomicratioset"