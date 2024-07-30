#' From a Sorensen-Dice threshold dissimilarity matrix, generate an object of class "hclust"
#'
#' @param x an object of class "dist" with the Sorensen-Dice equivalence threshold dissimilarities matrix
#' @param onTheFlyDev character, name of the graphical device where to immediately display the resulting
#'   diagram. The appropriate names depend on the operating system. Defaults to \code{NULL} and then
#'   nothing is displayed
#' @param method character, one of the admissible methods in function \code{hclust}. Defaults to "complete"
#' @param jobName character, main plot name, defaults to
#'   \code{paste("Equivalence cluster", onto, ontoLevel, method, sep = "_")}
#' @param ylab character, label of the vertical axis of the plot, defaults to "Sorensen equivalence threshold dissimilarity"
#' @param ... additional arguments to \code{hclust}
#'
#' @return An object of class \code{equivClustSorensen}, descending from class \code{hclust}
#'
#' @importFrom stringr str_pad
#' @importFrom stats hclust
#'
#' @examples
#' # Gene lists to analyse:
#' data("allOncoGeneLists")
#' 
#' # Obtaining ENTREZ identifiers for the gene universe of humans:
#' library(org.Hs.eg.db)
#' humanEntrezIDs <- keys(org.Hs.eg.db, keytype = "ENTREZID")
#' 
#' # First, compute the Sorensen-Dice threshold equivalence dissimilarity
#' # for ontology BP at level 4:
#' # # Very time consuming, it requires building all joint enrichment contingency tables
#' dOncBP4 <- sorenThreshold(allOncoGeneLists, onto = "BP", GOLevel = 4,
#'                           geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
#' # Better (much faster), using the previously tabulated contingency tables:
#' data("allTabsBP.4")
#' dOncBP4 <- sorenThreshold(allTabsBP.4)
#' clust.threshold <- hclustThreshold(dOncBP4)
#' plot(clust.threshold, main = "AllOnco genelists, BP ontology at level 4",
#'      ylab = "Sorensen equivalence threshold")
#' # With the same data, an UPGMA dendrogram:
#' clust.threshold <- hclustThreshold(dOncBP4, method = "average")
#' plot(clust.threshold, main = "AllOnco genelists, BP ontology at level 4",
#'      ylab = "Sorensen equivalence threshold")
#' @export
hclustThreshold <- function(x, onTheFlyDev = NULL, method = "complete",
                            jobName = paste("Equivalence cluster", method, sep = "_"),
                            ylab = "Sorensen equivalence threshold dissimilarity",
                            ...)
{
  subName <- paste0("Ontology ", attr(x, "onto"), " at level ", attr(x, "GOLevel"), sep = "")
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
