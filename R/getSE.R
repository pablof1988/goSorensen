#' Access to the estimated standard error of the sample Sorensen-Dice dissimilarity
#' in one or more equivalence test results
#'
#' Given objects representing the result(s) of one or more equivalence tests
#' (classes "equivSDhtest", "equivSDhtestList" or "allEquivSDtest", i.e., the
#' result of functions 'equivTestSorensen' and 'allEquivTestSorensen')
#' this function returns the estimated standard errors of the sample dissimilarities
#' in the tests.
#'
#' @param x an object of class "equivSDhtest" or "equivSDhtestList" or "allEquivSDtest".
#' @param onto character, a vector with one or more of "BP", "CC" or "MF", ontologies to access.
#' @param GOLevel numeric, a vector with one or more GO levels to access.
#' @param simplify logical, if TRUE the result is simplified, e.g., returning a vector instead
#' of a matrix.
#' @param listNames character(2), the names of a pair of gene lists.
#' @param ... Additional parameters.
#'
#' @ return A numeric value, the standard error of the Sorensen-Dice dissimilarity estimate
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
#' # (But results may vary according to GO updating)
#' getSE(waldman_atlas.BP.4)
#'
#' # All pairwise equivalence tests at level 4 of the BP ontology:
#' BP.4
#' class(BP.4)
#' # This may correspond to a call like:
#' # BP.4 <- equivTestSorensen(allOncoGeneLists,
#' #                           geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#' #                           onto = "BP", GOLevel = 4)
#' getSE(BP.4)
#' getSE(BP.4, simplify = FALSE)
#'
#' # Equivalence test iterated over all GO ontologies and levels 3 to 10:
#' cancerEquivSorensen
#' class(cancerEquivSorensen)
#' # This may correspond to code like:
#' # (By default, the tests are iterated over all GO ontologies and for levels 3 to 10)
#' # cancerEquivSorensen <- allEquivTestSorensen(allOncoGeneLists,
#' #                                             geneUniverse = humanEntrezIDs,
#' #                                             orgPackg = "org.Hs.eg.db")
#' # All standard errors of the Sorensen-Dice dissimilarity estimates:
#' getSE(cancerEquivSorensen)
#' getSE(cancerEquivSorensen, simplify = FALSE)
#' # Standard errors only for some GO ontologies, levels or pairs of gene lists:
#' getSE(cancerEquivSorensen, GOLevel = "level 6")
#' getSE(cancerEquivSorensen, GOLevel = "level 6", simplify = FALSE)
#' getSE(cancerEquivSorensen, GOLevel = "level 6", listNames = c("waldman", "sanger"))
#' getSE(cancerEquivSorensen, GOLevel = "level 6", onto = "BP")
#' getSE(cancerEquivSorensen, GOLevel = "level 6", onto = "BP", simplify = FALSE)
#' getSE(cancerEquivSorensen, GOLevel = "level 6", onto = "BP", listNames = c("waldman", "sanger"))
#' getSE(cancerEquivSorensen$BP$`level 4`)
#'

#'
#' @export
getSE <- function(x, ...) {
  UseMethod("getSE")
}

#' @describeIn getSE S3 method for class "equivSDhtest"
#' @export
getSE.equivSDhtest <- function(x, ...) {
  return(attr(x$estimate, "se"))
}

#' @describeIn getSE S3 method for class "equivSDhtestList"
#' @export
getSE.equivSDhtestList <- function(x, simplify = TRUE, ...) {
  result <- lapply(x, function(xi){
    resaux <- lapply(xi, getSE.equivSDhtest)
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
#' @describeIn getSE S3 method for class "AllEquivSDhtest"
#' @export
getSE.AllEquivSDhtest <- function(x, onto, GOLevel, listNames,
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
            return(attr(x[[ionto]][[ilev]][[ilist1]][[ilist2]]$estimate, "se"))
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
        return(attr(x[[ionto]][[ilev]][[listNames[1]]][[listNames[2]]]$estimate, "se"))
      }
    })
    names(resLev) <- GOLevel
    return(resLev)
  })
  names(result) <- onto
  return(result)
}
