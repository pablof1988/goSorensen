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
