#' Checks for validity data representing an enrichment contingency table generated from two gene lists
#'
#' @param x either an object of class "table", "matrix" or "numeric".
#'
#' @return boolean, TRUE if \code{x} nicely represents a 2x2 contingency table
#' interpretable as the cross-tabulation of the enriched GO items in two gene lists:
#' "Number of enriched items in list 1 (TRUE, FALSE)" x "Number of enriched items in
#' list 2 (TRUE, FALSE)". In this function, "nicely representing a 2x2 contingency table"
#' is interpreted in terms of computing the Sorensen-Dice dissimilarity and associated
#' statistics.
#' Otherwise the execution is interrupted.
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
#' tryCatch(nice2x2Table(conti), error = function(e){return(e)})
#' conti2 <- conti[1,seq.int(1, min(2,ncol(conti))), drop = FALSE]
#' conti2
#' tryCatch(nice2x2Table(conti2), error = function(e){return(e)})
#'
#' conti3 <- matrix(c(12, 210), ncol = 2, nrow = 1)
#' conti3
#' tryCatch(nice2x2Table(conti3), error = function(e){return(e)})
#'
#' conti4 <- c(32, 21, 81, 1439)
#' nice2x2Table(conti4)
#' conti4.mat <- matrix(conti4, nrow = 2)
#' conti4.mat
#' conti5 <- c(32, 21, 81)
#' nice2x2Table(conti5)
#'
#' conti6 <- c(-12, 21, 8)
#' tryCatch(nice2x2Table(conti6), error = function(e){return(e)})
#'
#' # Error: All enrichment frequencies are null: Inadequate for Sorensen-Dice computations
#' conti7 <- c(0, 0, 0, 32)
#' tryCatch(nice2x2Table(conti7), error = function(e){return(e)}) 


#'
#' @export
nice2x2Table <- function(x) {
  UseMethod("nice2x2Table")
}

#' @describeIn nice2x2Table S3 method for class "table"
#' @export
nice2x2Table.table <- function(x) {
  stopifnot("Not a 2x2 table" = dim(x) == c(2,2))
  stopifnot("Negative frequencies in a contingency table" = all(x >= 0))
  stopifnot("Zero frequencies: Inadequate for Sorensen-Dice computations" = sum(x[seq_len(3)]) > 0)
  return(TRUE)
}

#' @describeIn nice2x2Table S3 method for class "matrix"
#' @export
nice2x2Table.matrix <- function(x) {
  stopifnot("Not a 2x2 table" = dim(x) == c(2,2))
  stopifnot("Negative frequencies in a contingency table" = all(x >= 0))
  stopifnot("Zero frequencies: Inadequate for Sorensen-Dice computations" = sum(x[seq_len(3)]) > 0)
  return(TRUE)
}

#' @describeIn nice2x2Table S3 method for class "numeric"
#' @export
nice2x2Table.numeric <- function(x) {
  stopifnot("Not a 2x2 table" = dim(x) == c(2,2))
  stopifnot("Negative frequencies in a contingency table" = all(x >= 0))
  stopifnot("All enrichment frequencies are null: Inadequate for Sorensen-Dice computations" = sum(x[seq_len(3)]) > 0)
  return(TRUE)
}
