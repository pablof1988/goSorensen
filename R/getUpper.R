#' Access to the upper limit of the one-sided confidence intervals for the Sorensen-Dice
#' dissimilarity in one or more equivalence test results
#'
#' Given objects representing the result(s) of one or more equivalence tests
#' (classes "equivSDhtest", "equivSDhtestList" or "allEquivSDtest", i.e., the
#' result of functions 'equivTestSorensen' and 'allEquivTestSorensen')
#' this function returns the upper limits of the one-sided confidence intervals
#' [0, dU] for the Sorensen-Dice dissimilarity.
#'
#' @param x an object of class "equivSDhtest" or "equivSDhtestList" or "allEquivSDtest".
#' @param onto character, a vector with one or more of "BP", "CC" or "MF", ontologies to access.
#' @param GOLevel numeric or character, a vector with one or more GO levels to access.
#' See the details section and the examples.
#' @param simplify logical, if TRUE the result is simplified, e.g., returning a vector instead
#' of a matrix.
#' @param listNames character(2), the names of a pair of gene lists.
#' @param ... Additional parameters.
#'
#' @return A numeric value, the upper limit of the one-sided confidence interval for the Sorensen-Dice
#' dissimilarity.
#' @return When \code{x} is an object of class "equivSDhtest" (i.e., the result of a single
#' equivalence test), the returned value is a single numeric value, the upper limit of the
#' one-sided confidence interval for the Sorensen-Dice dissimilarity.
#' For an object of class "equivSDhtestList" (i.e. all pairwise tests for a
#' set of gene lists), if \code{simplify = TRUE} (the default), the resulting value is a vector
#' with the upper limit of the one-sided confidence intervals in all those tests, or the symmetric
#' matrix of all these values if \code{simplify = TRUE}. If \code{x} is an object of class
#' "allEquivSDtest" (i.e., the test iterated along GO ontologies and levels), the preceding result
#' is returned in the form of a list along the ontologies, levels and pairs of gene lists specified
#' by the arguments \code{onto, GOlevel} and \code{listNames} (or all present in \code{x} for
#' missing arguments).
#'
#' @details
#' Argument \code{GOLevel} can be of class "character" or "numeric". In the first case, the GO
#' levels must be specified like \code{"level 6"} or \code{c("level 4", "level 5", "level 6")}
#' In the second case ("numeric"), the GO levels must be specified like\code{6} or \code{4:6}.
#'
#' @examples
#' # Dataset 'allOncoGeneLists' contains the result of the equivalence test between gene lists
#' # 'waldman' and 'atlas', at level 4 of the BP ontology:
#' data(waldman_atlas.BP.4)
#' waldman_atlas.BP.4
#' class(waldman_atlas.BP.4)
#' # This may correspond to the result of code like:
#' # waldman_atlas.BP.4 <- equivTestSorensen(
#' #   allOncoGeneLists[["waldman"]], allOncoGeneLists[["atlas"]],
#' #   geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#' #   onto = "BP", GOLevel = 4, listNames = c("waldman", "atlas"))
#' # (But results may vary according to GO updating)
#' getUpper(waldman_atlas.BP.4)
#'
#' # All pairwise equivalence tests at level 4 of the BP ontology:
#' data(BP.4)
#' ?BP.4
#' class(BP.4)
#' # This may correspond to a call like:
#' # BP.4 <- equivTestSorensen(allOncoGeneLists,
#' #                           geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#' #                           onto = "BP", GOLevel = 4)
#' getUpper(BP.4)
#' getUpper(BP.4, simplify = FALSE)
#'
#' # Equivalence test iterated over all GO ontologies and levels 3 to 10:
#' data(cancerEquivSorensen)
#' ?cancerEquivSorensen
#' class(cancerEquivSorensen)
#' # This may correspond to code like:
#' # (By default, the tests are iterated over all GO ontologies and for levels 3 to 10)
#' # cancerEquivSorensen <- allEquivTestSorensen(allOncoGeneLists,
#' #                                             geneUniverse = humanEntrezIDs,
#' #                                             orgPackg = "org.Hs.eg.db")
#' # All upper confidence limits for the Sorensen-Dice dissimilarities:
#' getUpper(cancerEquivSorensen)
#' getUpper(cancerEquivSorensen, simplify = FALSE)
#'
#' # Upper confidence limits only for some GO ontologies, levels or pairs of gene lists:
#' getUpper(cancerEquivSorensen, GOLevel = "level 6")
#' getUpper(cancerEquivSorensen, GOLevel = 6)
#' getUpper(cancerEquivSorensen, GOLevel = 4:6)
#' getUpper(cancerEquivSorensen, GOLevel = "level 6", simplify = FALSE)
#' getUpper(cancerEquivSorensen, GOLevel = "level 6", listNames = c("waldman", "sanger"))
#' getUpper(cancerEquivSorensen, GOLevel = 4:6, onto = "BP")
#' getUpper(cancerEquivSorensen, GOLevel = 4:6, onto = "BP", simplify = FALSE)
#' getUpper(cancerEquivSorensen, GOLevel = "level 6", onto = "BP",
#'          listNames = c("waldman", "sanger"))
#' getUpper(cancerEquivSorensen$BP$`level 4`)

#'
#' @export
getUpper <- function(x, ...) {
  UseMethod("getUpper")
}

#' @describeIn getUpper S3 method for class "equivSDhtest"
#' @export
getUpper.equivSDhtest <- function(x, ...) {
  return(x$conf.int[2])
}

#' @describeIn getUpper S3 method for class "equivSDhtestList"
#' @export
getUpper.equivSDhtestList <- function(x, simplify = TRUE, ...) {
  result <- lapply(x, function(xi){
    resaux <- lapply(xi, getUpper.equivSDhtest)
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
  }else{
    return(result)
  }
}
#'
#' @describeIn getUpper S3 method for class "AllEquivSDhtest"
#' @export
getUpper.AllEquivSDhtest <- function(x, onto, GOLevel, listNames,
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
    stopifnot("'listNames' must be a 'character' of length 2" =
                is.character(listNames) && (length(listNames) == 2))
  }
  result <- lapply(onto, function(ionto) {
    resLev <- lapply(GOLevel, function(ilev) {
      if (allLists) {
        namsList1 <- names(x[[ionto]][[ilev]])
        namsMat <- c(names(x[[ionto]][[ilev]][[1]])[1], namsList1)
        resList1 <- sapply(namsList1, function(ilist1) {
          namsList2 <- names(x[[ionto]][[ilev]][[ilist1]])
          resList2 <- vapply(namsList2, function(ilist2) {
            return(x[[ionto]][[ilev]][[ilist1]][[ilist2]]$conf.int[2])
          }, FUN.VALUE = 0.0)
          names(resList2) <- namsList2
          return(resList2)
        })
        names(resList1) <- namsList1
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
      }else {
        return(x[[ionto]][[ilev]][[listNames[1]]][[listNames[2]]]$conf.int[2])
      }
    })
    names(resLev) <- GOLevel
    return(resLev)
  })
  names(result) <- onto
  return(result)
}
