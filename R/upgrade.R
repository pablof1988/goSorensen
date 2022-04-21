#' Update the result of a Sorensen-Dice equivalence test.
#'
#' Recompute the test (or tests) from an object of class "equivSDhtest", "equivSDhtestList" or
#' "AllEquivSDhtest" (i.e.,the output of functions "equivTestSorensen" or "allEquivTestSorensen").
#' Using the same table or tables of enrichment frequencies in 'x', obtain again the result of the
#' equivalence test for new values of any of the parameters \code{d0} or \code{conf.level} or
#' \code{boot} or \code{nboot} or \code{check.table}.
#'
#' @param x an object of class "equivSDhtest", "equivSDhtestList" or "AllEquivSDhtest"
#' @param ... any valid parameters for function "equivTestSorensen" for its interface "table".
#' The test(s) will be recomputed according to these parameters, or according to default values
#' for non-specified parameters
#'
#' @examples
#' # Result of the equivalence test between gene lists 'waldman' and 'atlas', in dataset
#' # 'allOncoGeneLists', at level 4 of the BP ontology:
#' waldman_atlas.BP.4
#' class(waldman_atlas.BP.4)
#' # This may correspond to the result of code like:
#' # waldman_atlas.BP.4 <- equivTestSorensen(
#' #   allOncoGeneLists[["waldman"]], allOncoGeneLists[["atlas"]],
#' #   geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#' #   onto = "BP", GOLevel = 4, listNames = c("waldman", "atlas"))
#' upgrade(waldman_atlas.BP.4, d0 = 1/(1 + 10/9)) # d0 = 0.4737
#' upgrade(waldman_atlas.BP.4, d0 = 1/(1 + 2*1.25)) # d0 = 0.2857
#' upgrade(waldman_atlas.BP.4, d0 = 1/(1 + 2*1.25), conf.level = 0.99)
#'
#' # All pairwise equivalence tests at level 4 of the BP ontology
#' BP.4
#' class(BP.4)
#' # This may correspond to a call like:
#' # BP.4 <- equivTestSorensen(allOncoGeneLists,
#' #                           geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#' #                           onto = "BP", GOLevel = 4)
#' upgrade(BP.4, d0 = 1/(1 + 2*1.25)) # d0 = 0.2857
#'
#' cancerEquivSorensen
#' class(cancerEquivSorensen)
#' upgrade(cancerEquivSorensen, d0 = 1/(1 + 2*1.25)) # d0 = 0.2857
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
