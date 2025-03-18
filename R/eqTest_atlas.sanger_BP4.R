#' Example of the output produced by the function \code{equivTestSorensen}. It contains the equivalence test for comparing two lists at level 4 of ontology BP.
#'
#' @description
#' The output of an equivalence test to detect biological similarity between the lists atlas and sanger from \code{allOncoGeneLists}, based on the normal asymptotic distribution.
#' 
#' @details
#' The parameters considered to execute this test are: irrelevance limit \code{d0 = 0.4444} and confidence level \code{conf.level = 0.95}.
#' 
#' Consider this object only as an illustrative example, which is valid exclusively for the lists atlas and sanger from the data \code{\link{allOncoGeneLists}} contained in this package. Note that gene lists, GO terms, and Bioconductor may change over time. The current version of these results were generated with Bioconductor version 3.20.
#' 
#' @format An exclusive object from `goSorensen` of the class "equivSDhtest" 
#' @usage data(eqTest_atlas.sanger_BP4)
"eqTest_atlas.sanger_BP4"