#' Remove all NULL or unrepresentable as a dendrogram "equivClustSorensen" elements in
#' an object of class "equivClustSorensenList"

#' @param x An object of class "equivClustSorensenList" descending from "iterEquivClust" which itself
#' descends from class "list". See the details section.

#' @return
#' An object of class "equivClustSorensenList".
#' @details
#' "equivClustSorensenList" objects are lists whose components are one or more of BP, CC or MF,
#' the GO ontologies. Each of these elements is itself a list whose elements correspond to GO levels.
#' Finally, the elements of these lists are objects of class "equivClustSorensen", descending
#' from "equivClust" which itself descends from "hclust".

#' @export
pruneClusts <- function(x) {
  result <- lapply(x, function(thisOntoClusts) {
    thisOntoClusts <- thisOntoClusts[vapply(thisOntoClusts,
                                            function(clust){
                                              !is.null(clust) &
                                                (length(attr(clust, "distMat")) >= 3)
                                            },
                                            FUN.VALUE = FALSE)]
  })
  class(result) <- class(x)
  return(result)
}
