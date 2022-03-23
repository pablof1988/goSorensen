
#' Creates a 2x2 enrichment contingency table from two gene lists
#'
#' @param x an object of class "character" (or coerzable to "character") representing a vector of gene identifiers.
#' @param y an object of class "character" (or coerzable to "character") representing a vector of gene identifiers.
#' @param listNames a character(2) object containing the names of the gene lists
#' which will originate the cross-tabulated enrichment frequencies.
#' @param check.table boolean. The resulting table must be checked. Defaults to TRUE.
#' @param geneUniverse character vector containing all genes from where geneLists have been extracted.
#' @param orgPackg A string with the name of the annotation package.
#' @param onto string describing the ontology. Either "BP", "MF" or "CC".
#' @param GOLevel An integer, the GO ontology level.
#' @param restricted Boolean variable to decide how tabulation of GOIDs is performed. Defaults to FALSE.
#' Unrestricted tabulation crosses _all_ GO Terms located at the level indicated by `GOLev`
#' with the two GOIDs lists. Restricted tabulation crosses only terms from the selected GO level
#' that are _common to ancestor terms of either list_.
#' That is, if one term in the selected GO level is not an ancestor of at least one of the gene list
#' most specific GO terms it is excluded from the GO Level's terms because it is impossible that it
#' appears as being enriched.
#' @param pAdjustMeth string describing the adjust method, either "BH", "BY" or "Bonf", defaults to 'BH'.
#' @param pvalCutoff A numeric value. Defaults to 0.05.
#' @param qvalCutoff A numeric value. Defaults to 0.01.
#' @import clusterProfiler goProfiles devtools GO.db org.Hs.eg.db
#'
#' @return an object of class "table" representing a 2x2 contingency table
#' interpretable as the cross-tabulation of the enriched GO items in two gene lists:
#' "Number of enriched items in list 1 (TRUE, FALSE)" x "Number of enriched items in
#' list 2 (TRUE, FALSE)".
#'
#' @examples
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
                                     check.table = TRUE, geneUniverse, orgPackg, onto, GOLevel,
                                     restricted = FALSE,
                                     pAdjustMeth = "BH", pvalCutoff = 0.01, qvalCutoff = 0.05)
{
  buildEnrichTable.character(as.character(x), as.character(y), listNames, check.table,
                             geneUniverse, orgPackg, onto, GOLevel,
                             restricted,
                             pAdjustMeth, pvalCutoff, qvalCutoff)
}

#' @describeIn buildEnrichTable S3 method for class "character"
#' @export
buildEnrichTable.character <- function(x, y, listNames = c("gene.list1", "gene.list2"),
                                       check.table = TRUE, geneUniverse, orgPackg, onto, GOLevel,
                                       restricted =FALSE,
                                       pAdjustMeth="BH", pvalCutoff=0.01, qvalCutoff=0.05) {
  tab <- crossTabGOIDs4GeneLists (genelist1 = x, genelist2 = y,
                                  geneUniverse, orgPackg, onto, GOLevel,
                                  restricted,
                                  pAdjustMeth, pvalCutoff, qvalCutoff)
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
