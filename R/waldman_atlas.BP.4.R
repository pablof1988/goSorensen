#' An example of "equivSDhtest" object resulting from a call to function 'equivSorensenTest'
#'
#' The Sorensen-Dice equivalence test between the gene lists "waldman" and "atlas" taken from dataset
#' \code{\link{allOncoGeneLists}} which is automatically charged with this package.
#' To perform the test, the information in these gene
#' lists was summarized by means of contingency tables of mutual GO item enrichment, for all GO items
#' at level 4 of the BP ontology. The tests were performed for an equivalence
#' limit d0 = 0.4444 and a confidence level conf.int = 0.95.
#' Based on a non up-to-date version of these gene lists, take just as an illustrative example.
#'
#' @format An object of class "equivSDhtest" inheriting from class "list".
#' @source \url{http://www.bushmanlab.org/links/genelists}
"waldman_atlas.BP.4"
