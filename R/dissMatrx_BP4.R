#' Example of the output produced by the function \code{sorenThreshold}. It contains the dissimilarity matrix at GO level 4, for the ontology BP.
#'
#' @description
#' This object contains the matrix of dissimilarities between the 7 lists from \code{allOncoGeneLists}, computed based on the irrelevance threshold that makes them equivalent at GO level 4, for the ontology BP. 
#' 
#' @details
#' Equivalence tests were computed based on the normal distribution (\code{boot = TRUE} by default) and using a confidence level \code{conf.level = 0.95}.
#' 
#' Consider this object only as an illustrative example, which is valid exclusively for the data \code{\link{allOncoGeneLists}} contained in this package. Note that gene lists, GO terms, and Bioconductor may change over time. The current version of these results were generated with Bioconductor version 3.20.
#' 
#' @format An object of class "dist" 
#' @usage data("dissMatrx_BP4")
"dissMatrx_BP4"