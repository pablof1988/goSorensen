
#' Creates a 2x2 enrichment contingency table from two gene lists
#'
#' @param x an object of class "character" (or coerzable to "character") representing a vector of gene identifiers.
#' @param y an object of class "character" (or coerzable to "character") representing a vector of gene identifiers.
#' @param listNames a character(2) with the gene lists names originating the cross-tabulated enrichment frequencies.
#' @param check.table boolean. The resulting table must be checked. Defaults to TRUE.
#' @param geneUniverse character vector containing all genes from where geneLists have been extracted.
#' @param orgPackg A string with the name of the annotation package.
#' @param onto string describing the ontology. Either "BP", "MF" or "CC".
#' @param GOLevel An integer, the GO ontology level.
#' @param restricted Boolean variable to decide how tabulation of GOIDs is performed. Defaults to FALSE.
#' See the details section.
#' @param pAdjustMeth string describing the adjust method, either "BH", "BY" or "Bonf", defaults to 'BH'.
#' @param pvalCutoff A numeric value. Defaults to 0.05.
#' @param qvalCutoff A numeric value. Defaults to 0.01.
#' @param ... Additional parameters for internal use (not used for the moment)
#'
#' @import clusterProfiler goProfiles devtools GO.db org.Hs.eg.db
#'
#' @return an object of class "table" representing a 2x2 contingency table
#' interpretable as the cross-tabulation of the enriched GO items in two gene lists:
#' "Number of enriched items in list 1 (TRUE, FALSE)" x "Number of enriched items in
#' list 2 (TRUE, FALSE)".
#'
#' @details
#' Unrestricted tabulation crosses _all_ GO Terms located at the level indicated by `GOLev`
#' with the two GOIDs lists. Restricted tabulation crosses only terms from the selected GO level
#' that are _common to ancestor terms of either list_.
#' That is, if one term in the selected GO level is not an ancestor of at least one of the gene list
#' most specific GO terms it is excluded from the GO Level's terms because it is impossible that it
#' appears as being enriched.

#' @examples
#' # Gene universe:
#' data(humanEntrezIDs)
#' # Gene lists to be explored for enrichment:
#' ?allOncoGeneLists
#' # Table of mutual GO node enrichment between gene lists Vogelstein and sanger,
#' # for ontology MF at GO level 6 (only first 50 genes, to improve speed).
#' vog.VS.sang <- buildEnrichTable(allOncoGeneLists[["Vogelstein"]][1:50],
#'                                 allOncoGeneLists[["sanger"]][1:50],
#'                                 geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#'                                 onto = "MF", GOLevel = 6, listNames = c("Vogelstein", "sanger"))
#' vog.VS.sang
#' # This is an inadequate table for Sorensen-Dice computations:
#' equivTestSorensen(vog.VS.sang)
#' # This sometimes happens, due too small gene lists or due to poor incidence
#' # of enrichment.
#' #
#' # In fact, the complete gene lists generate a much interesting contingency table:
#' # vog.VS.sang <- buildEnrichTable(allOncoGeneLists[["Vogelstein"]],
#' #                                 allOncoGeneLists[["sanger"]],
#' #                                 geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#' #                                 onto = "MF", GOLevel = 6, listNames = c("Vogelstein", "sanger"))
#' # vog.VS.sang
#' # equivTestSorensen(vog.VS.sang)

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
                                     pAdjustMeth = "BH", pvalCutoff = 0.01, qvalCutoff = 0.05, ...)
{
  buildEnrichTable.character(as.character(x), as.character(y), listNames, check.table,
                             geneUniverse, orgPackg, onto, GOLevel,
                             restricted,
                             pAdjustMeth, pvalCutoff, qvalCutoff, ...)
}

#' @describeIn buildEnrichTable S3 method for class "character"
#' @export
buildEnrichTable.character <- function(x, y, listNames = c("gene.list1", "gene.list2"),
                                       check.table = TRUE, geneUniverse, orgPackg, onto, GOLevel,
                                       restricted =FALSE,
                                       pAdjustMeth="BH", pvalCutoff=0.01, qvalCutoff=0.05, ...) {
  tab <- crossTabGOIDs4GeneLists (genelist1 = x, genelist2 = y,
                                  geneUniverse, orgPackg, onto, GOLevel,
                                  restricted,
                                  pAdjustMeth, pvalCutoff, qvalCutoff, ...)
  if (!all(dim(tab) == 2)) {
    tab <- completeTable(tab)
  }
  tab <- tab[c(2,1),c(2,1)]
  if (check.table){
    if (!nice2x2Table.table(tab)) {
      print(tab)
      stop("Inadequate GO terms enrichment contingency table")}
    if (sum(tab[1:3]) == 0) {
      warning("Zero enrichment frequencies: Inadequate table for Sorensen-Dice computations")
    }

  }
  dimnames(tab) = list(c(TRUE, FALSE), c(TRUE, FALSE))
  names(dimnames(tab)) <- paste0("Enriched in ", listNames)
  return(tab)
}
