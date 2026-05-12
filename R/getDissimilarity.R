#' Access to the estimated Sorensen-Dice dissimilarity in one or more e
#' quivalence test results
#'
#' Given objects representing the result(s) of one or more equivalence tests
#' (classes "equivSDhtest", "equivSDhtestList" or "allEquivSDtest", i.e., the
#' result of functions 'equivTestSorensen' and 'allEquivTestSorensen')
#' this function returns the estimated dissimilarities in the tests.
#'
#' @param x an object of class "equivSDhtest" or "equivSDhtestList" or
#' "allEquivSDtest".
#' @param onto character, a vector with one or more of "BP", "CC" or "MF",
#' ontologies to access.
#' @param GOLevel numeric or character, a vector with one or more GO levels to
#' access.
#' See the details section and the examples.
#' @param simplify logical, if TRUE the result is simplified, e.g., returning a
#' vector instead
#' of a matrix.
#' @param listNames character(2), the names of a pair of gene lists.
#' @param ... Additional parameters.
#'
#' @return When \code{x} is an object of class "equivSDhtest" (i.e., the result
#' of a single equivalence test), the returned value is a single numeric value,
#' the Sorensen-Dice dissimilarity. For an object of class "equivSDhtestList"
#' (i.e. all pairwise tests for a set of gene lists), if
#' \code{simplify = TRUE} (the default), the resulting value is a vector
#' with the dissimilarities in all those tests, or the symmetric matrix of all
#' dissimilarities if \code{simplify = TRUE}. If \code{x} is an object of class
#' "allEquivSDtest" (i.e., the test iterated along GO ontologies and levels),
#' the preceding result is returned in the form of a list along the ontologies,
#' levels and pairs of gene lists specified by the arguments
#' \code{onto, GOlevel} and \code{listNames} (or all present in \code{x} for
#' missing arguments).
#'
#' @details
#' Argument \code{GOLevel} can be of class "character" or "numeric". In the
#' first case, the GO levels must be specified like \code{"level 6"} or
#' \code{c("level 4", "level 5", "level 6")} In the second case ("numeric"),
#' the GO levels must be specified like\code{6} or \code{seq.int(4,6)}.
#'
#' @examples
#' # Manually define a 2 x 2 enrichment contingency table
#' contTable <- as.table(matrix(c(127, 19, 159, 3018),
#'   nrow = 2,
#'   dimnames = list(
#'     "Enriched in List 1" = c(TRUE, FALSE),
#'     "Enriched in List 2" = c(TRUE, FALSE)
#'   )
#' ))
#' contTable
#'
#' # Compute the Sorensen equivalence test from the contingency table.
#' # The result is an object of class "equivSDhtest".
#' equivalenceTest <- equivTestSorensen(contTable)
#' equivalenceTest
#'
#' # Extract the Sorensen dissimilarity value from the equivalence test object
#' getDissimilarity(equivalenceTest)
#'
#' @export
getDissimilarity <- function(x, ...) {
  UseMethod("getDissimilarity")
}

#' @describeIn getDissimilarity S3 method for class "equivSDhtest"
#' @export
getDissimilarity.equivSDhtest <- function(x, ...) {
  return(x$estimate)
}

#' @describeIn getDissimilarity S3 method for class "equivSDhtestList"
#' @export
getDissimilarity.equivSDhtestList <- function(x, simplify = TRUE, ...) {
  result <- lapply(x, function(xi) {
    resaux <- lapply(xi, getDissimilarity.equivSDhtest)
    names(resaux) <- names(xi)
    return(resaux)
  })
  namsList1 <- names(result) <- names(x)
  result <- unlist(result)
  if (!simplify) {
    namsMat <- c(names(x[[1]])[1], namsList1)
    nrows <- length(namsList1) + 1
    resMat <- matrix(0.0, nrow = nrows, ncol = nrows)
    resMat[upper.tri(resMat)] <- result
    resMat[lower.tri(resMat)] <- t(resMat)[lower.tri(resMat)]
    rownames(resMat) <- namsMat
    colnames(resMat) <- namsMat
    return(resMat)
  } else {
    return(result)
  }
}
#'
#' @describeIn getDissimilarity S3 method for class "AllEquivSDhtest"
#' @export
getDissimilarity.AllEquivSDhtest <- function(x, onto, GOLevel, listNames,
                                             simplify = TRUE, ...) {
  if (missing(onto)) {
    onto <- names(x)
  }
  if (missing(GOLevel)) {
    GOLevel <- names(x[[1]])
  } else {
    if (is.numeric(GOLevel)) {
      GOLevel <- paste0("level ", GOLevel)
    }
  }
  allLists <- missing(listNames)
  if (!allLists) {
    stopifnot(
      "'listNames' must be a 'character' of length 2" =
        is.character(listNames) && (length(listNames) == 2)
    )
  }
  result <- lapply(onto, function(ionto) {
    resLev <- lapply(GOLevel, function(ilev) {
      if (allLists) {
        namsList1 <- names(x[[ionto]][[ilev]])
        namsMat <- c(names(x[[ionto]][[ilev]][[1]])[1], namsList1)
        resList1 <- sapply(namsList1, function(ilist1) {
          namsList2 <- names(x[[ionto]][[ilev]][[ilist1]])
          resList2 <- vapply(namsList2, function(ilist2) {
            return(x[[ionto]][[ilev]][[ilist1]][[ilist2]]$estimate)
          }, FUN.VALUE = 0.0)
          return(resList2)
        })
        resList1 <- unlist(resList1, recursive = FALSE)
        if (!simplify) {
          nrows <- length(namsList1) + 1
          matList1 <- matrix(0.0, nrow = nrows, ncol = nrows)
          matList1[upper.tri(matList1)] <- resList1
          matList1[lower.tri(matList1)] <- t(matList1)[lower.tri(matList1)]
          rownames(matList1) <- namsMat
          colnames(matList1) <- namsMat
          resList1 <- matList1
        }
        return(resList1)
      } else {
        return(x[[ionto]][[ilev]][[listNames[1]]][[listNames[2]]]$estimate)
      }
    })
    names(resLev) <- GOLevel
    return(resLev)
  })
  names(result) <- onto
  return(result)
}
