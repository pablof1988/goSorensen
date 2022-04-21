#' Access to the effective number of bootstrap replicates in one or more equivalence test
#' results (only in their bootstrap version)
#'
#' Given objects representing the result(s) of one or more equivalence tests
#' (classes "equivSDhtest", "equivSDhtestList" or "allEquivSDtest", i.e., the
#' result of functions 'equivTestSorensen' and 'allEquivTestSorensen' with the parameter
#' boot = TRUE, otherwise it returns a NA value), this function returns the number of
#' effective bootstrap replicates.
#'
#' @param x an object of class "equivSDhtest" or "equivSDhtestList" or "allEquivSDtest".
#' @param onto character, a vector with one or more of "BP", "CC" or "MF", ontologies to access.
#' @param GOLevel numeric, a vector with one or more GO levels to access.
#' @param simplify logical, if TRUE the result is simplified, e.g., returning a vector instead
#' of a matrix.
#' @param listNames character(2), the names of a pair of gene lists.
#' @param ... Additional parameters.
#'
#' @examples
#' # Dataset 'allOncoGeneLists' contains the result of the equivalence test between gene lists
#' # 'waldman' and 'atlas', at level 4 of the BP ontology:
#' waldman_atlas.BP.4
#' class(waldman_atlas.BP.4)
#' # This may correspond to the result of code like:
#' # waldman_atlas.BP.4 <- equivTestSorensen(
#' #   allOncoGeneLists[["waldman"]], allOncoGeneLists[["atlas"]],
#' #   geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#' #   onto = "BP", GOLevel = 4, listNames = c("waldman", "atlas"))
#' #
#' # (But results may vary according to GO updating)
#'
#' # Not a bootstrap test, first upgrade to a bootstrap test:
#' boot.waldman_atlas.BP.4 <- upgrade(waldman_atlas.BP.4, boot = TRUE)
#'
#' getNboot(waldman_atlas.BP.4)
#' getNboot(boot.waldman_atlas.BP.4)
#'
#' getDissimilarity(waldman_atlas.BP.4)
#' getSE(waldman_atlas.BP.4)
#' getTable(waldman_atlas.BP.4)
#' getUpper(waldman_atlas.BP.4)
#' getDissimilarity(boot.waldman_atlas.BP.4)
#' getSE(boot.waldman_atlas.BP.4)
#' getTable(boot.waldman_atlas.BP.4)
#' getUpper(boot.waldman_atlas.BP.4)
#'
#'
#' # All pairwise equivalence tests at level 4 of the BP ontology
#' BP.4
#' class(BP.4)
#' # This may correspond to a call like:
#' # BP.4 <- equivTestSorensen(allOncoGeneLists,
#' #                           geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#' #                           onto = "BP", GOLevel = 4)
#' boot.BP.4 <- upgrade(BP.4, boot = TRUE)
#' getNboot(BP.4)
#' getNboot(boot.BP.4)
#' getNboot(boot.BP.4, simplify = FALSE)
#'
#' # Equivalence test iterated over all GO ontologies and levels 3 to 10.
#' cancerEquivSorensen
#' class(cancerEquivSorensen)
#' # This may correspond to code like:
#' # (By default, the tests are iterated over all GO ontologies and for levels 3 to 10)
#' # cancerEquivSorensen <- allEquivTestSorensen(allOncoGeneLists,
#' #                                             geneUniverse = humanEntrezIDs,
#' #                                             orgPackg = "org.Hs.eg.db")
#' boot.cancerEquivSorensen <- upgrade(cancerEquivSorensen, boot = TRUE)
#'
#' getNboot(cancerEquivSorensen)
#' getNboot(boot.cancerEquivSorensen)
#' getNboot(boot.cancerEquivSorensen, simplify = FALSE)
#' getNboot(boot.cancerEquivSorensen, GOLevel = "level 6")
#' getNboot(boot.cancerEquivSorensen, GOLevel = "level 6", simplify = FALSE)
#' getNboot(boot.cancerEquivSorensen, GOLevel = "level 6", listNames = c("waldman", "sanger"))
#' getNboot(boot.cancerEquivSorensen, GOLevel = "level 6", onto = "BP")
#' getNboot(boot.cancerEquivSorensen, GOLevel = "level 6", onto = "BP", simplify = FALSE)
#' getNboot(boot.cancerEquivSorensen, GOLevel = "level 6", onto = "BP",
#'          listNames = c("waldman", "sanger"))
#' getNboot(boot.cancerEquivSorensen$BP$`level 4`)

#'
#' @export
getNboot <- function(x, ...) {
  UseMethod("getNboot")
}

#' @describeIn getNboot S3 method for class "equivSDhtest"
#' @export
getNboot.equivSDhtest <- function(x, ...) {
  result <- attr(x$meth, "nboot")
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
getNboot.AllEquivSDhtest <- function(x, onto, GOLevel, listNames,
                                     simplify = TRUE, ...) {
  if (missing(onto)) {
    onto <- names(x)
  }
  if (missing(GOLevel)) {
    GOLevel <- names(x[[1]])
  }
  allLists <- missing(listNames)
  result <- lapply(onto, function(ionto) {
    resLev <- lapply(GOLevel, function(ilev) {
      if (allLists) {
        namsList1 <- names(x[[ionto]][[ilev]])
        namsMat <- c(names(x[[ionto]][[ilev]][[1]])[1], namsList1)
        resList1 <- sapply(namsList1, function(ilist1) {
          namsList2 <- names(x[[ionto]][[ilev]][[ilist1]])
          resList2 <- sapply(namsList2, function(ilist2) {
            # return(attr(x[[ionto]][[ilev]][[ilist1]][[ilist2]]$method, "nboot"))
            return(getNboot(x[[ionto]][[ilev]][[ilist1]][[ilist2]]))
          })
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
        # return(attr(x[[ionto]][[ilev]][[listNames[1]]][[listNames[2]]]$method, "nboot"))
        return(getNboot(x[[ionto]][[ilev]][[listNames[1]]][[listNames[2]]]))
      }
    })
    names(resLev) <- GOLevel
    # if (simplify) {
    #   resLev <- unlist(resLev, recursive = FALSE)
    # }
    return(resLev)
  })
  names(result) <- onto
  # if (simplify) {
  #   result <- unlist(result, recursive = FALSE)
  # }
  return(result)
}
