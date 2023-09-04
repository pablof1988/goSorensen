#' Iterate \code{sorenThreshold} along the specified GO ontologies and GO levels
#'
#' @param x either an object of class "list" or an object of class "allTableList". In the first
#' case, each of idata:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAWElEQVR42mNgGPTAxsZmJsVqQApgmGw1yApwKcQiT7phRBuCzzCSDSHGMKINIeDNmWQlA2IigKJwIssQkHdINgxfmBBtGDEBS3KCxBc7pMQgMYE5c/AXPwAwSX4lV3pTWwAAAABJRU5ErkJggg==ts elements must be a "character" vector of gene identifiers. In the second case,
#' the object corresponds to all contingency tables of joint enrichment along one or more GO
#' ontologies and one or more GO levels.
#' @param geneUniverse character vector containing all genes from where geneLists have been extracted
#' @param orgPackg a string with the name of the annotation package
#' @param boot boolean. If TRUE, the confidence intervals and the test p-values are computed by means
#' of a bootstrap approach instead of the asymptotic normal approach. Defaults to FALSE.
#' @param nboot numeric, number of initially planned bootstrap replicates. Ignored if
#' \code{boot == FALSE}. Defaults to 10000.
#' @param boot.seed starting random seed for all bootstrap iterations. Defaults to 6551.
#'   see the details section
#' @param ontos "character", GO ontologies to analyse.
#' @param GOLevels "integer", GO levels to analyse inside each one of these GO ontologies.
#' @param trace Logical. If TRUE (default), the (usually very time consuming) process
#' is traced along the specified GO ontologies and levels.
#' @param alpha simultaneous nominal significance level for the equivalence tests to be repeteadly performed,
#'   defaults to 0.05
#' @param precis numerical precision in the iterative search of the equivalence threshold dissimilarities,
#' @param ... extra parameters for function \code{buildEnrichTable}.
#'
#' @return
#' An object of class "distList". It is a list with as many components as GO ontologies have been
#' analysed. Each of these elements is itself a list with as many components as GO levels have been
#' analysed. Finally, the elements of these lists are objects of class "dist" with the Sorensen-Dice
#' equivalence threshold dissimilarity.
#'
#' @examples
#' # # This example is extremely time consuming, it scans two GO ontologies and three
#' # # GO levels inside them to perform the equivalence test.
#' # # Gene universe:
#' # data("humanEntrezIDs")
#' # # Gene lists to be explored for enrichment:
#' # data("allOncoGeneLists")
#' # allSorenThreshold(allOncoGeneLists,
#' #                   geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#' #
#' # Much faster:
#' # Object \code{allTabs} of class "allTableList" contains all the pairwise contingency tables of
#' # joint enrichment for the gene lists in \code{allOncoGeneLists}, obtained along all three GO
#' # ontologies and along GO levels 3 to 10:
#' data(allTabs)
#' dSors <- allSorenThreshold(allTabs, ontos = c("MF", "BP"), GOLevels = seq.int(4,6))
#' dSors$BP$`level 5`
#'
#' @export
allSorenThreshold <- function(x, ...) {
  UseMethod("allSorenThreshold")
}

#' @describeIn allSorenThreshold S3 method for class "list"
#' @export
allSorenThreshold.list <- function(x, geneUniverse, orgPackg,
                                   boot = FALSE, nboot = 10000, boot.seed = 6551,
                                   ontos = c("BP", "CC", "MF"), GOLevels = seq.int(3,10),
                                   trace = TRUE, alpha = 0.05, precis = 0.001,
                                   ...)
{
  tabs <- allBuildEnrichTable(x, ontos = ontos, GOLevels = GOLevels, trace = trace, ...)
  return(allSorenThreshold.allTableList(tabs, ontos = ontos, GOLevels = GOLevels, trace = trace))
}

#' @describeIn allSorenThreshold S3 method for class "allTableList"
#' @export
allSorenThreshold.allTableList <- function(x,
                                           boot = FALSE, nboot = 10000, boot.seed = 6551,
                                           ontos, GOLevels,
                                           trace = TRUE, alpha = 0.05, precis = 0.001,
                                           ...)
{
  if (missing(ontos)) {
    ontos <- names(x)
  }
  missGOLevels <- missing(GOLevels)
  allOntos <- lapply(ontos, function(ontoName) {
    onto <- x[[ontoName]]
    if (is.null(onto)) {
      return(NULL)
    }
    if (trace) {
      cat("\nOntology ", ontoName, "\n")
    }
    if (missGOLevels) {
      levNames <- names(onto)
    } else {
      levNames <- paste0("level ", GOLevels)
    }
    thisOnto <- lapply(levNames, function(levName) {
      GOlevel <- onto[[levName]]
      if (is.null((GOlevel))) {
        return(NULL)
      }
      if (trace) {
        cat("\n ", levName, "\n")
      }
      result <- sorenThreshold.tableList(GOlevel, onto = ontoName, GOLevel = levName,
                                  boot = boot, nboot = nboot, boot.seed = boot.seed,
                                  trace = trace, alpha = alpha, precis = precis, ...)
      return(result)
    })
    names(thisOnto) <- levNames
    return(thisOnto)
  })
  names(allOntos) <- ontos
  class(allOntos) <- c("distList", "list")
  return(allOntos)
}
