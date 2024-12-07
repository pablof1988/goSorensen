completeTable <- function (x, listNames)
{
  fullTab <- matrix(0, ncol = 2, nrow = 2,
                    dimnames = list(c(FALSE, TRUE), c(FALSE, TRUE)))
  fullTab[rownames(x), colnames(x)] <- x[rownames(x), colnames(x)]
  if (!missing(listNames)) {
    names(dimnames(fullTab)) <- paste0("Enriched in ", listNames)
  }
  return(fullTab)
}
