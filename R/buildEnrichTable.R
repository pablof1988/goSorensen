#' Creates a 2x2 enrichment contingency table from two gene lists
#'
#' @param x an object of class "character" (or coerzable to "character") representing a vector of gene identifiers.
#' @param y an object of class "character" (or coerzable to "character") representing a vector of gene identifiers.
#' @param listNames a character(2) object containing the names of the gene lists
#' which will originate the cross-tabulated enrichment frequencies.
#' @param check.table boolean. The resulting table must be checked
#' @param ... extra parameters for function \code{crossTabGOIDs4GeneLists} in package \code{equivStandardTest}.
#'
#' @return an object of class "table" representing a 2x2 contingency table
#' interpretable as the cross-tabulation of the enriched GO items in two gene lists:
#' "Number of enriched items in list 1 (TRUE, FALSE)" x "Number of enriched items in
#' list 2 (TRUE, FALSE)".
#'
#' @details
#'
#' @examples
#'
#' # Gene universe:
#' data(humanEntrezIDs)
#' # Gene lists to be explored for enrichment:
#' ?allOncoGeneLists
#' ?pbtGeneLists
#'
#' vog.VS.Wald <- buildEnrichTable(allOncoGeneLists[["Vogelstein"]], allOncoGeneLists[["waldman"]],
#'                                 listNames = c("Vogelstein","Vogelstein"),
#'                                 geneUniverse = humanEntrezIDs,
#'                                 orgPackg = "org.Hs.eg.db",
#'                                 onto = "BP", GOLevel = 4)
#' vog.VS.Wald
#' # Incomplete 2x2 table due to zero frequencies (no annotated items for list IRITD5):
#' buildEnrichTable(pbtGeneLists[["IRITD3"]], pbtGeneLists[["IRITD5"]],
#'                  geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db", onto = "MF", GOLevel = 7)
#' buildEnrichTable(pbtGeneLists[["IRITD3"]], pbtGeneLists[["IRITD5"]],
#'                  check.table = FALSE,
#'                  geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db", onto = "MF", GOLevel = 7)
#'
#' @export
buildEnrichTable <- function(x, ...) {
  UseMethod("buildEnrichTable")
}

#' @describeIn buildEnrichTable S3 default method
#' @export
buildEnrichTable.default <- function(x, y,
                                     listNames = c("gene.list1", "gene.list2"),
                                     check.table = TRUE, ...)
{
  buildEnrichTable.character(as.character(x), as.character(y), listNames, check.table, ...)
}

#' @describeIn buildEnrichTable S3 method for class "character"
#' @export
buildEnrichTable.character <- function(x, y, listNames = c("gene.list1", "gene.list2"),
                                       check.table = TRUE, ...) {
  tab <- crossTabGOIDs4GeneLists (genelist1 = x, genelist2 = y, ...)
  if (!all(dim(tab) == 2)) {
    tab <- completeTable(tab)
  }
  tab <- tab[c(2,1),c(2,1)]
  if (check.table){
    if (!nice2x2Table.table(tab)) {
      print(tab)
      stop("Inadequate GO items enrichment contingency table")}
  }
  dimnames(tab) = list(c(TRUE, FALSE), c(TRUE, FALSE))
  names(dimnames(tab)) <- paste0("Enriched in ", listNames)
  return(tab)
}
