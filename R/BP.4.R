#' An example of "equivSDhtestList" object resulting from a call to 'equivSorensenTest'
#'
#' All pairwise Sorensen-Dice equivalence tests between the gene lists in data \code{\link{allOncoGeneLists}} which
#' is automatically charged with this package.
#' To perform the tests, the information in these data was summarized as 2x2 contingency tables of GO items
#' enrichment, at level 4 of the BP ontology, and the tests were performed for an equivalence
#' limit d0 = 0.4444 and a confidence level conf.int = 0.95.
#'
#' @format An object of class "equivSDhtestList" inheriting from class "list". It is a list of class "equivSDhtest"
#' objects.
#' @source \url{http://www.bushmanlab.org/links/genelists}
"BP.4"
