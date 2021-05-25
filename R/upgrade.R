#' Update the result of an Sorensen-Dice equivalence test.
#'
#' Change d0 or conf.level in an object of class "equivSDhtest", "equivSDhtestList" or "AllEquivSDhtest"
#' (i.e.,the output of functions "equivTestSorensen" or "allEquivTestSorensen")
#'
#' @param x an object of class "equivSDhtest", "equivSDhtestList" or "AllEquivSDhtest"
#' @param ... additional parameters for functions "equivTestSorensen" or "allEquivTestSorensen",
#' for the moment, these are 'd0', the equivalence limit for the Sorensen-Dice dissimilarity,
#' and 'conf.level' the confidence level of the one-sided confidence interval upon which the
#' test is based.
#'
#' @examples
#' # Result of the equivalence test between gene lists 'waldman' and 'atlas', in dataset
#''allOncoGeneLists', at level 4 of the BP ontology:
#' waldman_atlas.BP.4
#' class(waldman_atlas.BP.4)
#' # This may correspond to the result of code like:
#' # waldman_atlas.BP.4 <- equivTestSorensen(allOncoGeneLists[["waldman"]], allOncoGeneLists[["atlas"]],
#' #                                         geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#' #                                         onto = "BP", GOLevel = 4, listNames = c("waldman", "atlas"))
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
  return(equivTestSorensen(tab, ...))
}

#' @describeIn upgrade S3 method for class "equivSDhtestList"
#' @export
upgrade.equivSDhtestList <- function(x, ...) {
  result <- lapply(x, function(xi){
    resaux <- lapply(xi, function(testij) {
      tab <- getTable.equivSDhtest(testij)
      return(equivTestSorensen(tab, ...))
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
          tab <- x[[ionto]][[ilev]][[ilist1]][[ilist2]]$enrichTab
          return(equivTestSorensen(tab, ...))
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
