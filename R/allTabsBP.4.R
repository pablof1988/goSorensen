#' An example of "tableList" object resulting from a call to 'buildEnrichTable'
#'
#' The result of generating all contingency tables of mutual enrichment, in a pairwise fashion,
#' between the gene lists in data \code{\link{allOncoGeneLists}}. The information in these data
#' was summarized as 2x2 contingency tables of GO items enrichment, at level 4 of the BP
#' ontology. These results are based on gene lists which are non automatically updated, take
#' them just as an illustrative example because the gene lists, and the GO, may change along
#' time.
#'
#' @format An object of class "tableList" inheriting from class "list". It is a list of class
#' "table" objects.
#' @usage data(allTabsBP.4)
"allTabsBP.4"
