#' Iterate \code{sorenThreshold} along the specified GO ontologies and GO levels
#'
#' @param x either an object of class "list" or an object of class "allTableList". In the first
#' case, each of its elements must be a "character" vector of gene identifiers (e.g., ENTREZ). In the second case,
#' the object corresponds to all contingency tables of joint enrichment along one or more GO
#' ontologies and one or more GO levels.
#' @param orgPackg A string with the name of the genomic annotation package corresponding to a specific species to be analyzed, which must be previously installed and activated. For more details see \href{../doc/README.html}{README File}.
#' @param geneUniverse character vector containing the universe of genes from where gene lists have been extracted. This vector must be obtained from the annotation package declared in orgPackg. For more details see \href{../doc/README.html}{README File}.
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
#' # # This example is highly time-consuming. It scans two GO ontologies and three
#' # # GO levels inside them to perform the equivalence test.
#' 
#' # Obtaining ENTREZ identifiers for the gene universe of humans:
#' # library(org.Hs.eg.db)
#' # humanEntrezIDs <- keys(org.Hs.eg.db, keytype = "ENTREZID")
#' 
#' # # Gene lists to be explored for enrichment:
#' # data("allOncoGeneLists")
#' # allSorenThreshold(allOncoGeneLists,
#' #                   geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#' #
#' # Much faster:
#' # Object allContTabs of class "allTableList" contains all the pairwise contingency tables of
#' # joint enrichment for the gene lists in \code{allOncoGeneLists}, obtained along all three GO
#' # ontologies and along GO levels 3 to 10:
#' data(allContTabs)
#' dSors <- allSorenThreshold(allContTabs, ontos = c("MF", "BP"), GOLevels = seq.int(4,6))
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
  tabs <- allBuildEnrichTable(x, geneUniverse = geneUniverse, orgPackg = orgPackg, 
                              boot = boot, nboot = nboot, boot.seed = boot.seed, 
                              ontos = ontos, GOLevels = GOLevels, trace = trace, 
                              alpha = alpha, precis = precis, ...)
  return(allSorenThreshold.allTableList(tabs, geneUniverse = geneUniverse, orgPackg = orgPackg,
                                        boot = boot, nboot = nboot, boot.seed = boot.seed,
                                        ontos = ontos, GOLevels = GOLevels, trace = trace,
                                        alpha = alpha, precis = precis, ...))
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
