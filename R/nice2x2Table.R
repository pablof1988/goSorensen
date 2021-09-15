#' Checks for validity data representing an enrichment contingency table generated from two gene lists
#'
#' @param x either an object of class "table", "matrix" or "numeric".
# @param n total number of enriched GO nodes. Only required (sometimes) on the "numeric" interface,
# see the details section.
#'
#' @return boolean, TRUE if \code{x} nicely represents a 2x2 contingency table
#' interpretable as the cross-tabulation of the enriched GO items in two gene lists:
#' "Number of enriched items in list 1 (TRUE, FALSE)" x "Number of enriched items in
#' list 2 (TRUE, FALSE)". In this function, "nicely representing a 2x2 contingency table"
#' is interpreted in terms of computing the Sorensen-Dice dissimilarity and associated
#' statistics. In these computations, the double negatives n00 are ignored.
#' The result is FALSE otherwise and warnings are issued.
#'
#' @details
#' In the "table" and "matrix" interfaces, the input parameter \code{x} must correspond
#' to a two-dimensional array:
#' \deqn{
#'  \tabular{rr}{
#'   n_{11} \tab n_{01} \cr
#'   n_{10} \tab n_{00},
#'  }
#' }{}
#'
#'\tabular{rr}{
#' n_11 \tab n_01 \cr
#' n_10 \tab n_00,
#'}
#' These values are interpreted (always in this order) as n11: number of GO items enriched in both lists,
#' n01: GO items enriched in the second list but not in the first one, n10: items not enriched in the second
#' list but enriched in the first one and double negatives, n00.
#' The double negatives n00 are ignored in many computations concerning the Sorensen-Dice index.
#'
#' In the "numeric" interface, the input \code{x} must correspond to a numeric of length
#' 3 or more, in the same order as before.
#'
#' @examples
#' conti <- as.table(matrix(c(27, 36, 12, 501, 43, 15, 0, 0, 0), nrow = 3, ncol = 3,
#'                          dimnames = list(c("a1","a2","a3"),
#'                                          c("b1", "b2","b3"))))
#' nice2x2Table(conti)
#' conti2 <- conti[1,1:min(2,ncol(conti)), drop = FALSE]
#' conti2
#' nice2x2Table(conti2)
#'
#' conti3 <- matrix(c(12, 210), ncol = 2, nrow = 1)
#' conti3
#' nice2x2Table(conti3)
#'
#' conti4 <- c(32, 21, 81, 1439)
#' nice2x2Table(conti4)
#' conti4.mat <- matrix(conti4, nrow = 2)
#' conti4.mat
#' conti5 <- c(32, 21, 81)
#' nice2x2Table(conti5)
#'


#'
#' @export
nice2x2Table <- function(x, ...) {
  UseMethod("nice2x2Table")
}

#' @describeIn nice2x2Table S3 method for class "table"
#' @export
nice2x2Table.table <- function(x) {
  if (!all(dim(x) == 2)) {
    message("Not a 2x2 contingency table")
    return(FALSE)
  }
  x123 <- x[1:3]
  if (sum(x123) == 0) {
    message("Inadequate frequencies for Sorensen-Dice computations")
    return(FALSE)
  }
  if (any(x123 < 0)) {
    message("Negative frequencies in a contingency table")
    return(FALSE)
  }
  return(TRUE)
}

#' @describeIn nice2x2Table S3 method for class "matrix"
#' @export
nice2x2Table.matrix <- function(x) {
  if (!all(dim(x) == 2)) {
    message("Not a 2x2 contingency table")
    return(FALSE)
  }
  x123 <- x[1:3]
  if (sum(x123) == 0) {
    message("Inadequate frequencies for Sorensen-Dice computations")
    return(FALSE)
  }
  if (any(x123 < 0)) {
    message("Negative frequencies in a contingency table")
    return(FALSE)
  }
  return(TRUE)
}

#' @describeIn nice2x2Table S3 method for class "numeric"
#' @export
nice2x2Table.numeric <- function(x) {
  if (length(x) < 3) {
    message("A numeric of almost length 3 is required to codify enrichment frequencies")
    return(FALSE)
  }
  x123 <- x[1:3]
  if (sum(x123) == 0) {
    message("Inadequate frequencies for Sorensen-Dice computations")
    return(FALSE)
  }
  if (any(x123 < 0)) {
    message("Negative frequencies in a contingency table")
    return(FALSE)
  }
  return(TRUE)
}
