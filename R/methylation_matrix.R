#' Matrix with simulated beta values for CpG islands
#' for testing purposes
#'
#' The matrix contains 20 rows (CpG islands) and 5 columns
#' (5 samples). The beta values are uniformly distributed
#' between 0 and 0.1 for all samples and CpG islands. 
#' Sample #1 values for CpG_4, _5, _6, _7, _8 are >> 0.1
#'
#' @usage data("methylation_matrix")
#' @return A matrix object.
#' @examples
#' data("methylation_matrix")
#' dim(methylation_matrix)
#'
"methylation_matrix"