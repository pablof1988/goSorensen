#' Iterate \code{equivTestSorensen} along the specified GO ontologies and GO levels
#'
#' @param x either an object of class "list" or an object of class "allTableList". In the first
#' case, each of its elements must be a "character" vector of gene identifiers (e.g., ENTREZ).
#' @param d0 equivalence threshold for the Sorensen-Dice dissimilarity, d.
#' The null hypothesis states that d >= d0, i.e., inequivalence between the compared
#' gene lists and the alternative that d < d0, i.e., equivalence or dissimilarity
#' irrelevance (up to a level d0).
#' @param conf.level confidence level of the one-sided confidence interval, a value between 0 and 1.
#' @param boot boolean. If TRUE, the confidence interval and the test p-value are computed by means
#' of a bootstrap approach instead of the asymptotic normal approach. Defaults to FALSE.
#' @param nboot numeric, number of initially planned bootstrap replicates. Ignored if
#' \code{boot == FALSE}. Defaults to 10000.
#' @param check.table Boolean. If TRUE (default), argument \code{x} is checked to adequately
#' represent a 2x2 contingency table (or an aggregate of them) or gene lists producing a correct
#' table. This checking is performed by means of function \code{nice2x2Table}.
#' @param ontos "character", GO ontologies to analyse. Defaults to \code{c("BP", "CC", "MF")}.
#' @param GOLevels "integer", GO levels to analyse inside each one of the GO ontologies.
#' @param trace Logical. If TRUE (default), the (usually very time consuming) process of function
#' \code{allEquivTestSorensen} is traced along the specified GO ontologies and levels.
#' @param ... extra parameters for function \code{buildEnrichTable}.
#'
#' @return
#' An object of class "AllEquivSDhtest". It is a list with as many components as GO ontologies have been analysed.
#' Each of these elements is itself a list with as many components as GO levels have been analized.
#' Finally, the elements of these lists are objects as generated by \code{equivTestSorensen.list},
#' i.e., objects of class "equivSDhtestList" containing pairwise comparisons between gene lists.
#'
#' @examples
#' # Gene lists to be explored for enrichment:
#' data(allOncoGeneLists)
#' 
#' # Obtaining ENTREZ identifiers for the gene universe of humans:
#' library(org.Hs.eg.db)
#' humanEntrezIDs <- keys(org.Hs.eg.db, keytype = "ENTREZID")
#' 
#' # This example is highly time-consuming. It scans two GO ontologies and three
#' # GO levels inside them to perform the equivalence test.
#' # allEquivTestSorensen(allOncoGeneLists,
#' #                      geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#' #                      ontos = c("MF", "BP"), GOLevels = seq.int(4,6))
#' # When the "ontos" and "GOLevels" arguments are not supplied, the function computes 
#' # by default every possible contingency table between the lists being compared for 
#' # the three ontologies (BP, CC, MF) and GO levels from 3 to 10. 
#' #
#' # Much faster:
#' # Object \code{allContTabs} of class "allTableList" contains all the pairwise contingency tables of
#' # joint enrichment for the gene lists in \code{allOncoGeneLists}, obtained along all three GO
#' # ontologies and along GO levels 3 to 10:
#' data(allContTabs)
#' tests <- allEquivTestSorensen(allContTabs, ontos = c("MF", "BP"), GOLevels = seq.int(4,6))
#' tests$BP$`level 5`
#' getPvalue(tests)
#'
#' @export
allEquivTestSorensen <- function(x, ...) {
  UseMethod("allEquivTestSorensen")
}

#' @describeIn allEquivTestSorensen S3 method for class "list"
#' @export
allEquivTestSorensen.list <- function(x, d0 = 1 / (1 + 1.25), conf.level = 0.95,
                                      boot = FALSE, nboot = 10000, check.table = TRUE,
                                      ontos = c("BP", "CC", "MF"),
                                      GOLevels = seq.int(3,10),
                                      trace = TRUE,
                                      ...)
{
  allOntos <- lapply(ontos, function(onto) {
    if (trace) {
      cat("\nPerforming all equivalence tests for ontology ", onto, "\n")
    }
    thisOnto <- lapply(GOLevels, function(lev) {
      if (trace) {
        cat("\n Level ", lev, "\n")
      }
      result <- equivTestSorensen(x, d0 = d0, conf.level = conf.level,
                                  boot = boot, nboot = nboot,
                                  check.table = check.table,
                                  onto = onto, GOLevel = lev, ...)
      return(result)
    })
    names(thisOnto) <- paste("level", GOLevels)
    return(thisOnto)
  })
  names(allOntos) <- ontos
  class(allOntos) <- c("AllEquivSDhtest", "list")
  return(allOntos)
}

#' @describeIn allEquivTestSorensen S3 method for class "allTableList"
#' @export
allEquivTestSorensen.allTableList <- function(x, d0 = 1 / (1 + 1.25), conf.level = 0.95,
                                              boot = FALSE, nboot = 10000, check.table = TRUE,
                                              ontos, GOLevels, trace = TRUE,
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
      cat("\nPerforming all equivalence tests for ontology ", ontoName, "\n")
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
        cat("\n", levName, "\n")
      }
      result <- equivTestSorensen(GOlevel, d0 = d0, conf.level = conf.level,
                                  boot = boot, nboot = nboot,
                                  check.table = check.table, ...)
      return(result)
    })
    names(thisOnto) <- levNames
    return(thisOnto)
  })
  names(allOntos) <- ontos
  class(allOntos) <- c("AllEquivSDhtest", "list")
  return(allOntos)
}
