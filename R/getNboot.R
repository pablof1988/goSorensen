#' Access to the number of initially planned bootstrap replicates in one or more equivalence test
#' results (only in their bootstrap version)
#'
#' Given objects representing the result(s) of one or more equivalence tests
#' (classes "equivSDhtest", "equivSDhtestList" or "allEquivSDtest", i.e., the
#' result of functions 'equivTestSorensen' and 'allEquivTestSorensen' with the parameter
#' boot = TRUE), this function returns the number of initially planned bootstrap replicates
#' in these equivalence tests, which may be greater than the number of finally effective or
#' valid bootstrap replicates. See the details section for more information on this.
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
#' @return When \code{x} is an object of class "equivSDhtest" (i.e., the result of a single
#' equivalence test), the returned value is a single numeric value, the number of initially
#' planned bootstrap replicates, or \code{NA} if bootstrapping has not been performed.
#' For an object of class "equivSDhtestList" (i.e. all pairwise tests for a
#' set of gene lists), if \code{simplify = TRUE} (the default), the resulting value is a vector
#' with the number of initially bootstrap replicates in all those tests, or the symmetric matrix
#' of all these values if \code{simplify = TRUE}.
#' If \code{x} is an object of class "allEquivSDtest"
#' (i.e., the test iterated along GO ontologies and levels), the preceding result is returned in
#' the form of a list along the ontologies, levels and pairs of gene lists specified by the arguments
#' \code{onto, GOlevel} and \code{listNames} (or all present in \code{x} for missing arguments).
#'
#' @details
#' In the bootstrap version of the equivalence test, resampling is performed generating new
#' bootstrap contingency tables from a multinomial distribution based on the "real", observed,
#' frequencies of mutual enrichment.
#' In some bootstrap iterations, the generated contingency table of mutual enrichment
#' may have very low frequencies of enrichment, which makes it unable for Sorensen-Dice
#' computations.
#' Then, the number of effective bootstrap resamples may be lower than those initially planned.
#' To get the number of effective bootstrap resamples use function \code{getEffNboot}.
#'
#' Argument \code{GOLevel} can be of class "character" or "numeric". In the first case, the GO
#' levels must be specified like \code{"level 6"} or \code{c("level 4", "level 5", "level 6")}
#' In the second case ("numeric"), the GO levels must be specified like\code{6} or \code{seq.int(4,6)}.
#'
#' @seealso \code{\link{getEffNboot}}
#'
#' @examples
#' # Dataset 'allOncoGeneLists' contains the result of the equivalence test between gene lists
#' # 'sanger' and 'atlas', at level 4 of the BP ontology:
#' data(eqTest_atlas.sanger_BP4)
#' eqTest_atlas.sanger_BP4
#' class(eqTest_atlas.sanger_BP4)
#' # This may correspond to the result of code like:
#' # eqTest_atlas.sanger_BP4 <- equivTestSorensen(
#' #   allOncoGeneLists[["sanger"]], allOncoGeneLists[["atlas"]],
#' #   geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#' #   onto = "BP", GOLevel = 4, listNames = c("sanger", "atlas"))
#' #
#' # (But results may vary according to GO updating)
#'
#' # Not a bootstrap test, first upgrade to a bootstrap test:
#' boot.eqTest_atlas.sanger_BP4 <- upgrade(eqTest_atlas.sanger_BP4, boot = TRUE)
#'
#' getNboot(eqTest_atlas.sanger_BP4)
#' getNboot(boot.eqTest_atlas.sanger_BP4)
#'
#' # All pairwise equivalence tests at level 4 of the BP ontology
#' data(eqTest_all_BP4)
#' ?eqTest_all_BP4
#' class(eqTest_all_BP4)
#' # This may correspond to a call like:
#' # eqTest_all_BP4 <- equivTestSorensen(allOncoGeneLists,
#' #                           geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#' #                           onto = "BP", GOLevel = 4)
#' boot.eqTest_all_BP4 <- upgrade(eqTest_all_BP4, boot = TRUE)
#' getNboot(eqTest_all_BP4)
#' getNboot(boot.eqTest_all_BP4)
#' getNboot(boot.eqTest_all_BP4, simplify = FALSE)
#'
#' # Bootstrap equivalence test iterated over all GO ontologies and levels 3 to 10.
#' # data(allEqTests)
#' # ?allEqTests
#' # class(allEqTests)
#' # This may correspond to code like:
#' # (By default, the tests are iterated over all GO ontologies and for levels 3 to 10)
#' # allEqTests <- allEquivTestSorensen(allOncoGeneLists,
#' #                                             geneUniverse = humanEntrezIDs,
#' #                                             orgPackg = "org.Hs.eg.db",
#' #                                             boot = TRUE)
#' # boot.allEqTests <- upgrade(allEqTests, boot = TRUE)
#' # All numbers of bootstrap replicates:
#' # getNboot(boot.allEqTests)
#' # getNboot(boot.allEqTests, simplify = FALSE)
#'
#' # Number of bootstrap replicates for specific GO ontologies, levels or pairs of gene lists:
#' # getNboot(boot.allEqTests, GOLevel = "level 6")
#' # getNboot(boot.allEqTests, GOLevel = 6)
#' # getNboot(boot.allEqTests, GOLevel = seq.int(4,6))
#' # getNboot(boot.allEqTests, GOLevel = "level 6", simplify = FALSE)
#' # getNboot(boot.allEqTests, GOLevel = "level 6", listNames = c("atlas", "sanger"))
#' # getNboot(boot.allEqTests, GOLevel = seq.int(4,6), onto = "BP")
#' # getNboot(boot.allEqTests, GOLevel = seq.int(4,6), onto = "BP", simplify = FALSE)
#' # getNboot(boot.allEqTests, GOLevel = "level 6", onto = "BP",
#' #          listNames = c("waldman", "sanger"))
#' # getNboot(boot.allEqTests$BP$`level 4`)

#'
#' @export
getNboot <- function(x, ...) {
  UseMethod("getNboot")
}

#' @describeIn getNboot S3 method for class "equivSDhtest"
#' @export
getNboot.equivSDhtest <- function(x, ...) {
  result <- attr(x, "nboot")
  return(if (is.null(result)) NA else result)
}

#' @describeIn getNboot S3 method for class "equivSDhtestList"
#' @export
getNboot.equivSDhtestList <- function(x, simplify = TRUE, ...) {
  result <- lapply(x, function(xi){
    resaux <- lapply(xi, getNboot.equivSDhtest)
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
#' @describeIn getNboot S3 method for class "AllEquivSDhtest"
#' @export
getNboot.AllEquivSDhtest <- function(x, onto, GOLevel, listNames, simplify = TRUE, ...) {
  if (missing(onto)) {
    onto <- names(x)
  }
  stopifnot("GO ontology names must be one or more of 'BP', 'CC' or 'MF'" =
              onto %in% c("BP","CC","MF"))
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
            # return(attr(x[[ionto]][[ilev]][[ilist1]][[ilist2]]$method, "nboot"))
            return(getNboot(x[[ionto]][[ilev]][[ilist1]][[ilist2]]))
          }, FUN.VALUE = 0)
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
        return(getNboot(x[[ionto]][[ilev]][[listNames[1]]][[listNames[2]]]))
      }
    })
    names(resLev) <- GOLevel
    return(resLev)
  })
  names(result) <- onto
  return(result)
}
