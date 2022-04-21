#' Access to the contingency table of mutual enrichment of one or more equivalence test results
#'
#' Given objects representing the result(s) of one or more equivalence tests
#' (classes "equivSDhtest", "equivSDhtestList" or "allEquivSDtest", i.e., the
#' result of functions 'equivTestSorensen' and 'allEquivTestSorensen')
#' this function returns the contingency tables from which the tests were performed.
#'
#' @param x an object of class "equivSDhtest" or "equivSDhtestList" or "allEquivSDtest".
#' @param onto character, a vector with one or more of "BP", "CC" or "MF", ontologies to access.
#' @param GOLevel numeric, a vector with one or more GO levels to access.
#' @param listNames character(2), the names of a pair of gene lists.
#' @param ... Additional parameters.

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
#' # (But results may vary according to GO updating)
#' getTable(waldman_atlas.BP.4)
#'
#' # All pairwise equivalence tests at level 4 of the BP ontology
#' BP.4
#' class(BP.4)
#' # This may correspond to a call like:
#' # BP.4 <- equivTestSorensen(allOncoGeneLists,
#' #                           geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#' #                           onto = "BP", GOLevel = 4)
#' getTable(BP.4)
#'
#' # Equivalence test iterated over all GO ontologies and levels 3 to 10:
#' cancerEquivSorensen
#' class(cancerEquivSorensen)
#' # This may correspond to code like:
#' # cancerEquivSorensen <- allEquivTestSorensen(allOncoGeneLists,
#' #                                             geneUniverse = humanEntrezIDs,
#' #                                             orgPackg = "org.Hs.eg.db")
#' # (By default, the tests are iterated over all GO ontologies and for levels 3 to 10)
#' # All 2x2 contingecy tables of joint enrichment:
#' getTable(cancerEquivSorensen)
#' # Contingency tables only for some GO ontologies, levels or pairs of gene lists:
#' getTable(cancerEquivSorensen, GOLevel = "level 6")
#' getTable(cancerEquivSorensen, GOLevel = "level 6", listNames = c("waldman", "sanger"))
#' getTable(cancerEquivSorensen, GOLevel = "level 6", onto = "BP")
#' getTable(cancerEquivSorensen, GOLevel = "level 6", onto = "BP",
#'          listNames = c("waldman", "sanger"))
#'
#'
#' @export
getTable <- function(x, ...) {
  UseMethod("getTable")
}

#' @describeIn getTable S3 method for class "equivSDhtest"
#' @export
getTable.equivSDhtest <- function(x, ...) {
  return(x$enrichTab)
}

#' @describeIn getTable S3 method for class "equivSDhtestList"
#' @export
getTable.equivSDhtestList <- function(x, ...) {
  result <- lapply(x, function(xi){
    resaux <- lapply(xi, getTable.equivSDhtest)
    names(resaux) <- names(xi)
    return(resaux)
  })
  names(result) <- names(x)
  return(result)
  # return(lapply(x, getTable.equivSDhtest))
}

#' @describeIn getTable S3 method for class "AllEquivSDhtest"
#' @export
getTable.AllEquivSDhtest <- function(x,
                                     onto, GOLevel, listNames, ...) {
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
        resList1 <- lapply(namsList1, function(ilist1) {
          namsList2 <- names(x[[ionto]][[ilev]][[ilist1]])
          resList2 <- lapply(namsList2, function(ilist2) {
            return(x[[ionto]][[ilev]][[ilist1]][[ilist2]]$enrichTab)
          })
          names(resList2) <- namsList2
          return(resList2)
        })
        names(resList1) <- namsList1
        # if (simplify) {
        #   resList1 <- unlist(resList1, recursive = FALSE)
        # }
        return(resList1)
      }else {
        return(x[[ionto]][[ilev]][[listNames[1]]][[listNames[2]]]$enrichTab)
      }
    })
    names(resLev) <- GOLevel
    # if (simplify && allLists) {
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
