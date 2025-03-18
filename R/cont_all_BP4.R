#' Example of the output produced by the function \code{buildEnrichTable}. It contains the enrichment contingency tables for all the lists from \code{allOncoGeneLists} at level 4 of ontology BP.
#'
#' @description
#' Given 7 lists contained  in \code{allOncoGeneLists}, this object contains the 7(6)/2 = 21 possible enrichment contingency tables to compare all possible pairs of lists.
#' Each contingency 2x2 table contains the number of joint enriched GO terms (TRUE-TRUE); the number of GO terms enriched only in one list but not in the other one (FALSE-TRUE and TRUE-FALSE); and the number of GO terms not enriched in either of the two lists.
#' 
#' An important attribute of this object is \code{enriched}, which contains the enrichment matrix obtained using the function \code{\link{enrichedIn}}.  Actually, the contingency tables in this object are derived from cross-frequency tables created between pairs of lists, which are located as columns in this enrichment matrix.
#' 
#' @details
#' Consider this object only as an illustrative example, which is valid exclusively for the data \code{\link{allOncoGeneLists}} contained in this package. Note that gene lists, GO terms, and Bioconductor may change over time. The current version of these results were generated with Bioconductor version 3.20.
#' 
#' @format An exclusive object from \code{goSorensen} of the class "tableList" 
#' @usage data(cont_all_BP4)
"cont_all_BP4"