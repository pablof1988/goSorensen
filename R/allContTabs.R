#' Example of the output produced by the function \code{allBuildEnrichTable}. 
#'
#' @description
#' This object contains all the enrichment contingency tables to compare all possible pairs of lists from \code{\link{allOncoGeneLists}} across GO-Levels 3 to 10, and for the ontologies BP, CC, and MF.
#' 
#' @details
#' The attribute \code{enriched} is present in each element of this output, meaning that there is an enrichment matrix, similar to the one obtained with the function \code{enrichedIn}, for each ontology and GO-Level contained in this object.
#' 
#' Consider this object only as an illustrative example, which is valid exclusively for the data \code{\link{allOncoGeneLists}} contained in this package. Note that gene lists, GO terms, and Bioconductor may change over time. The current version of these results were generated with Bioconductor version 3.20.
#' 
#' @format An exclusive object from \code{goSorensen} of the class "allTableList" 
#' @usage data(allContTabs)
"allContTabs"