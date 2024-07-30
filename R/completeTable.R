#' Reformats and completes (if necessary) a 2x2 enrichment contingency table for its appropriate use in package goSorensen. It is internally used by function 'buildEnrichTable'..
#' 
#' @param x an object of class "table".
#' @param listNames a character(2) with the gene lists names.
#' enrichment frequencies.
#' @return a complete contingency table to use in package goSorensen.
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
