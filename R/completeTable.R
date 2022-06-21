#' completeTable
#' 
#' Reformats and completes (if necessary) an enrichment contingency table as is generated
#' by function 'crossTabGOIDs4GeneLists', in order to make it appropriate for its use in
#' package goSorensen.
#'
#' @param x an object of class "table", typically the output of function 'crossTabGOIDs4GeneLists'.
#' @param listNames a character(2) with the gene lists names originating the cross-tabulated
#' enrichment frequencies.
completeTable <- function(x, listNames) {
  fullTab <- matrix(0, ncol = 2, nrow = 2,  dimnames = list(c(FALSE, TRUE),
                                                            c(FALSE, TRUE)))
  if (ncol(x) == 1) {
    if (nrow(x) == 1) {
      fullTab[1,1] <- x[1,1]
    } else {
      fullTab[,colnames(x)] <- x[,1]
    }
  } else {
    if (nrow(x) == 1) {
      fullTab[rownames(x),] <- x[1,]
    }
  }
  if (!missing(listNames)) {
    names(dimnames(fullTab)) <- paste0("Enriched in ", listNames)
  }
  return(fullTab)
}
