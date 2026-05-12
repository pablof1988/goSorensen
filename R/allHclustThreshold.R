#' Iterate \code{hclustThreshold} along the specified GO ontologies and GO
#' levels
#'
#' @param x an object of class "distList".
#' @param ontos "character", GO ontologies to iterate. Defaults to the
#' ontologies in 'x'.
#' @param GOLevels "integer", GO levels to iterate inside each one of these GO
#' ontologies.
#' @param trace Logical. If TRUE (default), the process
#' is traced along the specified GO ontologies and levels.
#' @param ... extra parameters for function \code{hclustThreshold}.
#'
#' @return
#' An object of class "equivClustSorensenList" descending from "iterEquivClust"
#' which itself descends
#' from class "list".
#' It is a list with as many components as GO ontologies have been
#' specified. Each of these elements is itself a list with as many components as
#' GO levels have been
#' specified. Finally, the elements of these lists are objects of class
#' "equivClustSorensen", descending
#' from "equivClust" which itself descends from "hclust".
#'
#' @examples
#' # The following example requires calculating all dissimilarity matrices at
#' # the 3:10 levels and BP, CC, and MF ontologies for visualization purposes.
#' # Because this process is computationally intensive and can take a
#' # considerable amount of time, the example is not run automatically during
#' # the R CMD check.
#' \dontrun{
#' ## i) Obtaining ENTREZ identifiers for the gene universe of humans:
#' library(org.Hs.eg.db)
#' humanEntrezIDs <- keys(org.Hs.eg.db, keytype = "ENTREZID")
#'
#' ## ii) Gene lists to be explored for analysis:
#' data(allOncoGeneLists)
#'
#' # iii) Compute the thresholded dissimilarity matrices for the BP, CC, and MF
#' # ontologies across all GO levels, specifically from level 3 to level 10.
#' allDissMatrx <- allSorenThreshold(allOncoGeneLists,
#'   geneUniverse = humanEntrezIDs,
#'   orgPackg = "org.Hs.eg.db",
#'   ontos = c("BP", "CC", "MF"),
#'   GOLevels = 3:10,
#'   trace = FALSE
#' )
#' allDissMatrx
#' }
#'
#' # Since running this example may take several minutes, the result has been
#' # pre-computed and is accessible as the following:
#' data(allDissMatrx)
#' allDissMatrx
#' # This shortcut applies only to this example; for your own gene-list data,
#' # the computation must be performed explicitly.
#'
#' # The hclusts object stores all hierarchical clustering results for each
#' # ontology and GO level calculated from the allDissMatrx object. Any of these
#' # clustering results can be visualized. For example:
#' all.clust.threshold <- allHclustThreshold(allDissMatrx)
#' plot(all.clust.threshold$BP$`level 4`)
#' plot(all.clust.threshold$CC$`level 5`)
#' plot(all.clust.threshold$MF$`level 6`)
#'
#' @export
allHclustThreshold <- function(x,
                               ontos, GOLevels,
                               trace = TRUE, ...) {
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
      message("\nOntology ", ontoName, "\n")
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
        message("\n ", levName, "\n")
      }
      result <- hclustThreshold(GOlevel,
        onto = ontoName, GOLevel = levName,
        trace = trace, ...
      )
      return(result)
    })
    names(thisOnto) <- levNames
    return(thisOnto)
  })
  names(allOntos) <- ontos
  class(allOntos) <- c("equivClustSorensenList", "iterEquivClust", "list")
  return(allOntos)
}
