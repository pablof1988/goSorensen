#' An example of an object of class "allTableList" resulting from a call to 'buildEnrichTable'
#'
#' The result of generating all contingency tables of mutual enrichment, in a pairwise fashion,
#' between the gene lists in data \code{\link{allOncoGeneLists}} along GO levels 3 to 10 for all
#' three GO ontologies, BP, MF and CC. Object 'allTabs' is a list of length 3 with
#' one element for each GO ontology: allTabs$BP, allTabs$MF and allTabs$CC. Each one of these lists
#' is itself a list of length 11, e.g., allTabs$BP$`level 3`, allTabs$BP$`level 4`, etc. to
#' allTabs$BP$`level 10`. Finally, each one of these lists contains all contingency tables of mutual
#' enrichment, in a pairwise fashion, for the gene lists in data \code{\link{allOncoGeneLists}}.
#' These results are based on gene lists which are not automatically
#' updated, take them just as an illustrative example because the gene lists, and the GO, may change along
#' time. The present version of these data was generated under Bioconductor version 3.17.
#'
#' @format An object of class "allTableList" inheriting from class "list".
#'
#' @usage data(allTabs)
"allTabs"
