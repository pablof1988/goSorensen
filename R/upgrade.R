#' Update the result of a Sorensen-Dice equivalence test.
#'
#' Recompute the test (or tests) from an object of class "equivSDhtest", "equivSDhtestList" or
#' "AllEquivSDhtest" (i.e.,the output of functions "equivTestSorensen" or "allEquivTestSorensen").
#' Using the same table or tables of enrichment frequencies in 'x', obtain again the result of the
#' equivalence test for new values of any of the parameters \code{d0} or \code{conf.level} or
#' \code{boot} or \code{nboot} or \code{check.table}.
#'
#' @param x an object of class "equivSDhtest", "equivSDhtestList" or "AllEquivSDhtest".
#' @param ... any valid parameters for function "equivTestSorensen" for its interface "table",
#' to recompute the test(s) according to these parameters.
#'
#' @return An object of the same class than \code{x}.
#'
#' @examples
#' # Result of the equivalence test between gene lists 'sanger' and 'atlas', in dataset
#' # 'allOncoGeneLists', at level 4 of the BP ontology:
#' data(eqTest_atlas.sanger_BP4)
#' eqTest_atlas.sanger_BP4
#' class(eqTest_atlas.sanger_BP4)
#' # This may correspond to the result of code like:
#' # data(allOncoGeneLists)
#' # library(org.Hs.eg.db)
#' # humanEntrezIDs <- keys(org.Hs.eg.db, keytype = "ENTREZID")
#' # eqTest_atlas.sanger_BP4 <- equivTestSorensen(
#' #   allOncoGeneLists[["sanger"]], allOncoGeneLists[["atlas"]],
#' #   geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#' #   onto = "BP", GOLevel = 4, listNames = c("sanger", "atlas"))
#' upgrade(eqTest_atlas.sanger_BP4, d0 = 1/(1 + 10/9)) # d0 = 0.4737
#' upgrade(eqTest_atlas.sanger_BP4, d0 = 1/(1 + 2*1.25)) # d0 = 0.2857
#' upgrade(eqTest_atlas.sanger_BP4, d0 = 1/(1 + 2*1.25), conf.level = 0.99)
#'
#' # All pairwise equivalence tests at level 4 of the BP ontology
#' data(eqTest_all_BP4)
#' ?eqTest_all_BP4
#' class(eqTest_all_BP4)
#' # This may correspond to a call like:
#' # data(allOncoGeneLists)
#' # library(org.Hs.eg.db)
#' # humanEntrezIDs <- keys(org.Hs.eg.db, keytype = "ENTREZID")
#' # eqTest_all_BP4 <- equivTestSorensen(allOncoGeneLists,
#' #                           geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#' #                           onto = "BP", GOLevel = 4)
#' upgrade(eqTest_all_BP4, d0 = 1/(1 + 2*1.25)) # d0 = 0.2857
#'
#' data(allEqTests)
#' ?allEqTests
#' class(allEqTests)
#' upgrade(allEqTests, d0 = 1/(1 + 2*1.25)) # d0 = 0.2857
#'
#' @export
upgrade <- function(x, ...) {
  UseMethod("upgrade")
}

#' @describeIn upgrade S3 method for class "equivSDhtest"
#' @export
upgrade.equivSDhtest <- function(x, ...) {
  tab <- getTable.equivSDhtest(x)
  newPars <- list(...)
  if (is.null(newPars$d0)) {
    d0 <- x$null.value
  } else {
    d0 <- newPars$d0
  }
  if (is.null(newPars$conf.level)) {
    conf.level <- attr(x$conf.int, "conf.level")
  } else {
    conf.level <- newPars$conf.level
  }
  oldBoot <- !is.null(attr(x$meth, "nboot"))
  if (is.null(newPars$boot)) {
    boot <- oldBoot
  } else {
    boot <- newPars$boot
  }
  if (is.null(newPars$nboot)) {
    if (oldBoot) {
      nboot <- attr(x$meth, "nboot")
    } else {
      if (boot) {
        nboot <- 10000
      } else {
        nboot <- NULL
      }
    }
  } else {
    nboot <- newPars$nboot
  }
  if (is.null(newPars$check.table)) {
    check.table <- attr(x, "check.table")
  } else {
    check.table <- newPars$check.table
  }
  return(equivTestSorensen.table(tab, d0 = d0, conf.level = conf.level,
                                 boot = boot, nboot = nboot, check.table = check.table))
}

#' @describeIn upgrade S3 method for class "equivSDhtestList"
#' @export
upgrade.equivSDhtestList <- function(x, ...) {
  result <- lapply(x, function(xi){
    resaux <- lapply(xi, function(testij) {
      tab <- getTable.equivSDhtest(testij)
      return(equivTestSorensen.table(tab, ...))
    })
    names(resaux) <- names(xi)
    return(resaux)
  })
  names(result) <- names(x)
  class(result) <- c("equivSDhtestList", "list")
  return(result)
}

#' @describeIn upgrade S3 method for class "allEquivSDhtest"
#' @export
upgrade.AllEquivSDhtest <- function(x, ...) {
  onto <- names(x)
  GOLevel <- names(x[[1]])
  result <- lapply(onto, function(ionto) {
    resLev <- lapply(GOLevel, function(ilev) {
      namsList1 <- names(x[[ionto]][[ilev]])
      resList1 <- lapply(namsList1, function(ilist1) {
        namsList2 <- names(x[[ionto]][[ilev]][[ilist1]])
        resList2 <- lapply(namsList2, function(ilist2) {
          tab <- getTable.equivSDhtest(x[[ionto]][[ilev]][[ilist1]][[ilist2]])
          return(equivTestSorensen.table(tab, ...))
        })
        names(resList2) <- namsList2
        return(resList2)
      })
      names(resList1) <- namsList1
      class(resList1) <- c("equivSDhtestList", "list")
      return(resList1)
    })
    names(resLev) <- GOLevel
    return(resLev)
  })
  names(result) <- onto
  class(result) <- c("AllEquivSDhtest", "list")
  return(result)
}
