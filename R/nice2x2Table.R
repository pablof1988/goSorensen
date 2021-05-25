#' Checks for validity data representing an enrichment contingency table and tries to put it
#' in an adequate format, if needed
#'
#' @param x either an object of class "table", "matrix" or "numeric".
#' @param listNames a character(2) object containing the names of the gene lists
#' which have originated the cross-tabulated enrichment frequencies. Ignore it to
#' let unchanged the original names (if any) in argument \code{x}.
#' @param n total number of cases (GO nodes) in the contingency table. Only required
#' (sometimes) on the "numeric" interface, see the details section.
#'
#' @return either an object of the same class than the input \code{x}
#' (i.e., "table", "matrix" or "numeric") nicely representing a 2x2 contingency table
#' interpretable as the cross-tabulation of the enriched GO items in two gene lists:
#' "Number of enriched items in list 1 (FALSE, TRUE)" x "Number of enriched items in
#' list 2 (FALSE, TRUE)".
#'
#' @details
#' In the "table" and "matrix" interfaces, the input parameter \code{x} must correspond
#' to a two-dimensional array. It is trimmed if one or both dimensions exceed 2 (and a
#' warning is issued) and, conversely, it is filled with zeros in order to complete
#' a 2x2 table, if required. This 2x2 table is finally returned.
#' In the "numeric" interface, the input \code{x} must correspond to a numeric of length
#' 3 or 4. In the first case (length 3), it is interpreted as the frequencies of:
#' "enriched in the first gene list and not in the second gene list",
#' "not enriched in the first list and enriched in the second" and
#' "enriched in both lists", in this order. Then it is completed to be a numeric of length
#' 4, adding at its first position the frequency of "not enriched in both lists", as
#' \code{n - sum(x)}. Otherwise the parameter \code{n} is ignored. This method always returns
#' a length 4 numeric.
#'
#' @examples
#' conti <- as.table(matrix(c(501, 27, 36, 12, 43, 15, 0, 0, 0),
#' nrow = 3, ncol = 3,
#' dimnames = list(c("a1","a2","a3"),
#'                 c("b1", "b2","b3"))))
#' nice2x2Table(conti)
#' nice2x2Table(conti, listNames = c("a gene list","another"))
#' conti2 <- conti[1,1:min(2,ncol(conti)), drop = FALSE]
#' conti2
#' nice2x2Table(conti2)
#' nice2x2Table(conti2, listNames = c("a gene list","another"))
#'
#' conti3 <- matrix(c(210, 12), ncol = 2, nrow = 1)
#' conti3
#' nice2x2Table(conti3)
#'
#' conti4 <- c(1439, 32, 21, 81)
#' nice2x2Table(conti4)
#' conti4.mat <- matrix(conti4, nrow = 2)
#' conti4.mat
#' conti5 <- c(32, 21, 81)
#' nice2x2Table(conti5, n = 1573)
#' nice2x2Table(conti5, n = 1000)
#' try(nice2x2Table(conti5, n = 10))
#'
#' # Building enrichment contingency tables from scratch, using package "equivStandardTest"
#' library(equivStandardTest)
#' # Gene universe:
#' data(humanEntrezIDs)
#' # Gene lists to be explored for enrichment:
#' data(kidneyGeneLists)
#' # Incomplete 2x2 table due to zero frequencies (no annotated items for list IRITD5):
#' IRITD3VS5.MF7 <- crossTabGOIDs4GeneLists(kidneyGeneLists[["IRITD3"]], kidneyGeneLists[["IRITD5"]],
#'   geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db", onto = "MF", GOLevel = 7)
#' IRITD3VS5.MF7
#' # Complete as a 2x2 table:
#' IRITD3VS5.MF7 <- nice2x2Table(IRITD3VS5.MF7, c("IRITD3", "IRITD5"))
#' IRITD3VS5.MF7


#'
#' @export
nice2x2Table <- function(x, ...) {
  UseMethod("nice2x2Table")
}

#' @describeIn nice2x2Table S3 method for class "table"
#' @export
nice2x2Table.table <- function(x, listNames) {
  if (!all(dim(x) == 2)) {
    if (any(dim(x) > 2)) {
      x <- x[1:min(2,dim(x)[1]), 1:min(2,dim(x)[2]), drop = FALSE]
      warning("One or more table dimensions are greater than 2, it has been trimmed")
    }
    if (any(dim(x) < 2)) {
      x <- as.table(completeTable(x))
    }
  }
  if (any(x < 0)) {
    stop("Negative frequencies in a contingency table")
  }
  if (!missing(listNames)) {
    dimnames(x) = list(c(FALSE, TRUE), c(FALSE, TRUE))
    names(dimnames(x)) <- paste0("Enriched in ", listNames)
  }
  return(x)
}

#' @describeIn nice2x2Table S3 method for class "matrix"
#' @export
nice2x2Table.matrix <- function(x, listNames) {
  if ((ncol(x) != 2) || (nrow(x) != 2)) {
    if (any(dim(x) > 2)) {
      x <- x[1:min(2,nrow(x)), 1:min(2,ncol(x)), drop = FALSE]
      warning("One or more table dimensions are greater than 2, it has been trimmed")
    }
    if (any(dim(x) < 2)) {
      x <- completeTable(x)
    }
  }
  if (any(x < 0)) {
    stop("Negative frequencies in a contingency table")
  }
  if (!missing(listNames)) {
    dimnames(x) = list(c(FALSE, TRUE), c(FALSE, TRUE))
    names(dimnames(x)) <- paste0("Enriched in ", listNames)
  }
  return(x)
}

#' @describeIn nice2x2Table S3 method for class "numeric"
#' @export
nice2x2Table.numeric <- function(x, n) {
  if (length(x) < 3) {
    stop("The minimum length of a vector representing a 2x2 enrichment contingency table must be 3")
  }
  if (length(x) == 3) {
    tab <- c(n - sum(x), x)
    if (!is.null(names(tab))) {
      names(tab)[1] <- "n00"
    } else {
      names(tab) <- c("n00", "n10", "n01", "n11")
    }
  } else {
    tab <- x[1:4]
    if (is.null(names(tab))) {
      names(tab) <- c("n00", "n10", "n01", "n11")
    }
  }
  if (any(tab < 0)) {
    stop("Negative frequencies in a contingency table")
  }
  return(tab)
}
