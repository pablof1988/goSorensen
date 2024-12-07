#' Iterate \code{hclustThreshold} along the specified GO ontologies and GO levels
#'
#' @param x an object of class "distList".
#' @param ontos "character", GO ontologies to iterate. Defaults to the ontologies in 'x'.
#' @param GOLevels "integer", GO levels to iterate inside each one of these GO ontologies.
#' @param trace Logical. If TRUE (default), the process
#' is traced along the specified GO ontologies and levels.
#' @param ... extra parameters for function \code{hclustThreshold}.
#'
#' @return
#' An object of class "equivClustSorensenList" descending from "iterEquivClust" which itself descends
#' from class "list".
#' It is a list with as many components as GO ontologies have been
#' specified. Each of these elements is itself a list with as many components as GO levels have been
#' specified. Finally, the elements of these lists are objects of class "equivClustSorensen", descending
#' from "equivClust" which itself descends from "hclust".
#'
#' @examples
#' # Object \code{allTabs} of class "allTableList" contains all the pairwise contingency tables of
#' # joint enrichment for the gene lists in \code{allOncoGeneLists}, obtained along all three GO
#' # ontologies and along GO levels 3 to 10:
#' data(allContTabs)
#' # Compute the Sorensen-Dice equivalence threshold dissimilarity (only for the MF and CC
#' # ontologies and from levels 4 to 6):
#' dists <- allSorenThreshold(allContTabs, ontos = c("MF", "CC"), GOLevels = seq.int(4,6))
#' hclusts <- allHclustThreshold(dists)
#' hclusts$MF$`level 6`
#' plot(hclusts$MF$`level 6`)
#'

#' @export
allHclustThreshold <- function(x,
                               ontos, GOLevels,
                               trace = TRUE, ...)
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
      result <- hclustThreshold(GOlevel, onto = ontoName, GOLevel = levName,
                                trace = trace, ...)
      return(result)
    })
    names(thisOnto) <- levNames
    return(thisOnto)
  })
  names(allOntos) <- ontos
  class(allOntos) <- c("equivClustSorensenList", "iterEquivClust", "list")
  return(allOntos)
}
