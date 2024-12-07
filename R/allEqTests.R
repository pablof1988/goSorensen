#' Example of the output produced by the function \code{allEquivTestSorensen} using the normal asymptotic distribution. 
#'
#' @description
#' This object contains all the outputs for the equivalence tests to compare all possible pairs of lists from \code{\link{allOncoGeneLists}} across GO-Levels 3 to 10, and for the ontologies BP, CC, and MF, using the normal asymptotic distribution.
#' 
#' @details
#' The parameters considered to execute these tests are: irrelevance limit \code{d0 = 0.4444} and confidence level \code{conf.level = 0.95}.
#' 
#' Consider this object only as an illustrative example, which is valid exclusively for the data \code{\link{allOncoGeneLists}} contained in this package. Note that gene lists, GO terms, and Bioconductor may change over time. The current version of these results were generated with Bioconductor version 3.20.
#' 
#' @format An exclusive object from \code{goSorensen} of the class "AllEquivSDhtest" 
#' @usage data(allEqTests)
"allEqTests"