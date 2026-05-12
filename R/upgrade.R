#' Update the result of a Sorensen-Dice equivalence test.
#'
#' Recompute the test (or tests) from an object of class "equivSDhtest",
#' "equivSDhtestList" or "AllEquivSDhtest" (i.e.,the output of functions
#' "equivTestSorensen" or "allEquivTestSorensen").
#' Using the same table or tables of enrichment frequencies in 'x', obtain again
#' the result of the equivalence test for new values of any of the parameters
#' \code{d0} or \code{conf.level} or \code{boot} or \code{nboot} or
#' \code{check.table}.
#'
#' @param x an object of class "equivSDhtest", "equivSDhtestList" or
#' "AllEquivSDhtest".
#' @param ... any valid parameters for function "equivTestSorensen" for its
#' interface "table", to recompute the test(s) according to these parameters.
#'
#' @return An object of the same class than \code{x}.
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
#' # Access the equivalence test results.
#' upgrade(equivalenceTest)
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
  return(equivTestSorensen.table(tab,
    d0 = d0, conf.level = conf.level,
    boot = boot, nboot = nboot,
    check.table = check.table
  ))
}

#' @describeIn upgrade S3 method for class "equivSDhtestList"
#' @export
upgrade.equivSDhtestList <- function(x, ...) {
  result <- lapply(x, function(xi) {
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
