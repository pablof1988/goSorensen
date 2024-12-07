#' Example of the output produced by the function \code{enrichedIn}. It contains all the GO terms enriched or not-enriched in the lists of \code{allOncoGeneLists}, ontology BP, GO-Level 4. 
#' 
#' @description
#' A matrix with columns representing the gene lists from \code{\link{allOncoGeneLists}}, and rows with GO terms in the BP ontology at GO-Level 4.
#' 
#' This matrix comprises logit values, with \code{TRUE} indicating that the associated GO term is enriched in the respective list, and \code{FALSE} indicating that the GO term is not enriched.
#' 
#' This matrix represents the output of the \code{\link{enrichedIn}} function with the argument \code{onlyEnriched = FALSE}. The rows of this matrix display all the GO terms involved in the BP ontology at GO-Level 4.
#' 
#' @details
#' The attribute \code{nTerms} indicates the total number of GO terms evaluated in the BP ontology, GO-Level 4. For this particular case, \code{nTerms} matches with the number of rows of the matrix
#' 
#' Please, consider this object as an illustrative example only, which is valid exclusively for the \code{\link{allOncoGeneLists}} data contained in this package. Please note that gene lists, GO terms and Bioconductor may change over time. The current version of these results was generated with Bioconductor version 3.20.
#' 
#' @format An object of class "matrix" "array" 
#' @usage data(fullEnrichedInBP4)
"fullEnrichedInBP4"