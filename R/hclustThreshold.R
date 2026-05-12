#' From a Sorensen-Dice threshold dissimilarity matrix, generate an object of c
#' lass "hclust"
#'
#' @param x an object of class "dist" with the Sorensen-Dice equivalence
#' threshold dissimilarities matrix
#' @param onTheFlyDev character, name of the graphical device where to
#' immediately display the resulting
#'   diagram. The appropriate names depend on the operating system. Defaults to
#'   \code{NULL} and then nothing is displayed
#' @param method character, one of the admissible methods in function
#' \code{hclust}. Defaults to "complete"
#' @param jobName character, main plot name, defaults to
#'   \code{paste("Equivalence cluster", onto, ontoLevel, method, sep = "_")}
#' @param ylab character, label of the vertical axis of the plot, defaults to
#' "Sorensen equivalence threshold dissimilarity"
#' @param ... additional arguments to \code{hclust}
#'
#' @return An object of class \code{equivClustSorensen}, descending from class
#' \code{hclust}
#'
#' @importFrom stringr str_pad
#' @importFrom stats hclust
#'
#' @examples
#'
#' # The following example requires the computation of the dissimilarity matrix
#' # for visualization purposes. Since this process is computationally intensive
#' # and may take a considerable amount of time, the example is not executed
#' # automatically during R CMD check
#'
#' \dontrun{
#' ## i) Obtaining ENTREZ identifiers for the gene universe of humans:
#' library(org.Hs.eg.db)
#' humanEntrezIDs <- keys(org.Hs.eg.db, keytype = "ENTREZID")
#'
#' ## ii) Gene lists to be explored for analysis:
#' data(allOncoGeneLists)
#'
#' # iii) Computing the threshold dissimilarity matrix directly from gene lists
#' # at GO level 4 and BP ontology:
#' dissMatrx_BP4 <- sorenThreshold(allOncoGeneLists,
#'   geneUniverse = humanEntrezIDs,
#'   orgPackg = "org.Hs.eg.db",
#'   onto = "BP",
#'   GOLevel = 4,
#'   trace = FALSE
#' )
#' dissMatrx_BP4
#' }
#'
#' # Since running this example may take several minutes, the result has been
#' # pre-computed and is accessible as the following:
#' data(dissMatrx_BP4)
#' dissMatrx_BP4
#' # This shortcut applies only to this example; for your own gene-list data,
#' # the computation must be performed explicitly.
#'
#' # Visualization of the hierarchical clustering results based on the
#' # dissimilarity matrix and irrelevance threshold values.
#' clust.threshold <- hclustThreshold(dissMatrx_BP4)
#' plot(clust.threshold)
#'
#' @export
hclustThreshold <- function(x, onTheFlyDev = NULL, method = "complete",
                            jobName = paste("Equivalence cluster",
                              method,
                              sep = "_"
                            ),
                            ylab = "Sorensen equivalence threshold
                            dissimilarity",
                            ...) {
  subName <- paste0("Ontology ", attr(x, "onto"), " at level ",
    attr(x, "GOLevel"),
    sep = ""
  )
  clust <- hclust(x, method = method)
  attr(clust, "jobName") <- jobName
  attr(clust, "ylab") <- ylab
  attr(clust, "sub") <- subName
  attr(clust, "distMat") <- x

  if (!is.null(onTheFlyDev) & (length(x) > 2)) {
    eval(call(onTheFlyDev, width = 20, height = 20))
    # dev.set(dev.list()[onTheFlyDev])
    plot(clust, hang = -1, main = jobName, sub = subName, ylab = ylab)
  }

  class(clust) <- c("equivClustSorensen", "equivClust", class(clust))
  return(clust)
}
