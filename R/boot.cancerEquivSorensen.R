#' An example of object of class "AllEquivSDhtest" resulting from a call to 'allEquivTestSorensen'
#'
#' The bootstrap Sorensen-Dice test performed on the cancer gene lists in data \code{\link{allOncoGeneLists}}
#' which is may be charged from this package. The test is iterated for all GO ontologies and for GO
#' levels 3 to 10. These results are not automatically updated for changes in these gene lists, take them just
#' as an illustrative example. The present version was obtained under Bioconductor 3.17.
#'
#' For each ontology and GO level, the result contains the result of all pairwise tests of equivalence between
#' the cancer gene lists.
#'
#' @format An object of class "AllEquivSDhtest" inheriting from class "list". Each one of its elements, named
#' BP, CC and MF respectively, corresponds to a GO ontology. It is itself a list of length 8 whose elements
#' are named as "Level 3" to "Level 10". For each combination of ontology and level, there is an object of
#' class "equivSDhtestList" codifying the result of all pairwise tests between these cancer gene lists.
#' @source \url{http://www.bushmanlab.org/links/genelists}
#' @usage data(boot.cancerEquivSorensen)
"boot.cancerEquivSorensen"
