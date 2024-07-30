#' An example of an object of class "equivSDhtest" resulting from a call to function 'equivSorensenTest'
#'
#' The Sorensen-Dice equivalence test between the gene lists "waldman" and "atlas" taken from dataset
#' \code{\link{allOncoGeneLists}} which may be charged from this package.
#' To perform the test, the information in these gene
#' lists was summarized by means of contingency tables of mutual GO term enrichment, for all GO terms
#' at level 4 of the BP ontology. The tests were performed for an equivalence
#' limit d0 = 0.4444 and a confidence level conf.int = 0.95.
#' Based on a version of these gene lists that may be non up-to-date, take just as an illustrative example.
#' The present result was obtained under Bioconductor 3.17.
#'
#' @format An object of class "equivSDhtest" inheriting from class "list".
#' @source \url{http://www.bushmanlab.org/links/genelists}
#' @usage data(waldman_atlas.BP.4)
"waldman_atlas.BP.4"
