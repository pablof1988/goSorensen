#' Standard error of the sample Sorensen-Dice dissimilarity, asymptotic approach
#'
#' @param x either an object of class "table", "matrix" or "numeric" representing a 2x2 contingency table,
#' or a "character" (a set of gene identifiers) or "list" object. See the details section for more information.
#' @param y an object of class "character" representing a vector of gene identifiers.
#' @param check.table Boolean. If TRUE (default), argument \code{x} is checked to adequately
#' represent a 2x2 contingency table. This checking is performed by means of function
#' \code{nice2x2Table}.
#' @param listNames character(2), names of both gene lists being compared (only for the character interface)
#' @param ... extra parameters for function \code{buildEnrichTable}.
#'
#' @return In the "table", "matrix", "numeric" and "character" interfaces, the value of the standard error of the
#' Sorensen-Dice dissimilarity estimate. In the "list" interface, the symmetric matrix of all pairwise standard error
#' dissimilarity estimates.
#'
#' @details
#' This function computes the standard error estimate of the sample Sorensen-Dice dissimilarity,
#' given a 2x2 arrangement of frequencies (either implemented as a "table", a "matrix"
#' or a "numeric" object):
#' \deqn{
#'  \tabular{rr}{
#'   n_{11} \tab n_{10} \cr
#'   n_{01} \tab n_{00},
#'  }
#' }{}
#'
#'\tabular{rr}{
#' n_11 \tab n_10 \cr
#' n_01 \tab n_00,
#'}
#'
#' The subindex '11' corresponds to those
#' GO items enriched in both lists, '01' to items enriched in the second list but not in the first one,
#' '10' to items enriched in the first list but not enriched in the second one and '00' corresponds
#' to those GO items non enriched in both gene lists, i.e., to the double negatives, a value which
#' is ignored in the computations.
#'
#' In the "numeric" interface, if \code{length(x) >= 3}, the values are interpreted
#' as
#' \eqn{(n_{11}, n_{01}, n_{10})}{%
#' (n_11, n_01, n_10)}, always in this order.
#'
#' If \code{x} is an object of class "character", it must represent a list of gene identifiers. Then the
#' standard error for the dissimilarity between lists \code{x} and \code{y} is computed, after summarizing
#' these gene lists as a 2x2 contingency table of joint enrichment.
#'
#' In the "list" interface, the argument must be a list of "character" vectors, each one
#' representing a gene list (character identifiers). Then, all pairwise standard errors
#' of the dissimilarity between these gene lists are computed.
#'
#' @seealso \code{\link{nice2x2Table}} for checking and reformatting data,
#' \code{\link{dSorensen}} for computing the Sorensen-Dice dissimilarity,
#' \code{\link{duppSorensen}} for the upper limit of a one-sided confidence interval
#' of the dissimilarity, \code{\link{equivTestSorensen}} for an equivalence test.
#'
#' @examples
#' # Gene lists 'atlas' and 'sanger' in 'Cangenes' dataset. Table of joint enrichment
#' # of GO items in ontology BP at level 3.
#' tab_atlas.sanger_BP3
#' dSorensen(tab_atlas.sanger_BP3)
#' seSorensen(tab_atlas.sanger_BP3)
#'
#' # Contingency table as a numeric vector:
#' seSorensen(c(56, 1, 30, 47))
#' seSorensen(c(56, 1, 30))
#'
#' # (These examples may be considerably time consuming due to many enrichment 
#' # tests to build the contingency tables of mutual enrichment)
#' # ?pbtGeneLists
#' # Standard error of the sample Sorensen-Dice dissimilarity, directly from
#' # two gene lists, from scratch:
#' # seSorensen(pbtGeneLists[[2]], pbtGeneLists[[4]],
#' #            listNames = names(pbtGeneLists)[c(2,4)],
#' #            onto = "CC", GOLevel = 5,
#' #            geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
#' # Essentially, the above code makes the same as:
#' # tab.IRITD3vsKT1 <- buildEnrichTable(pbtGeneLists[[2]], pbtGeneLists[[4]],
#' #                                     listNames = names(pbtGeneLists)[c(2,4)],
#' #                                     onto = "CC", GOLevel = 5,
#' #                                     geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
#' # tab.IRITD3vsKT1
#' # seSorensen(tab.IRITD3vsKT1)
#'
#' # All pairwise standard errors (quite time consuming):
#' # seSorensen(pbtGeneLists,
#' #            onto = "CC", GOLevel = 5,
#' #            geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
#'
#' @export
seSorensen <- function(x, ...) {
  UseMethod("seSorensen")
}

#' @describeIn seSorensen S3 method for class "table"
#' @export
seSorensen.table <- function(x,
                             listNames = NULL, check.table = TRUE, ...) {
  if (check.table){
    if (!nice2x2Table.table(x)) {
      print(x)
      stop("Inadequate table to compute the standard error")
    }
  }
  n <- sum(x[1:3])
  p11 <- x[1] / n
  p11plus <- 1 + p11
  result <- 2 * sqrt(p11 * (1 - p11) / (n - 1)) / (p11plus * p11plus)
  if (is.null(listNames)) {
    names(result) <- NULL
  } else {
    names(result) <- paste("Sorensen-Dice disimilarity standard error ", listNames[1], ",", listNames[2], sep = "")
  }
  return(result)
}

#' @describeIn seSorensen S3 method for class "matrix"
#' @export
seSorensen.matrix <- function(x, listNames = NULL, check.table = TRUE, ...) {
  if (check.table){
    if (!nice2x2Table.matrix(x)) {
      print(x)
      stop("Inadequate table to compute the standard error")
    }
  }
  n <- sum(x[1:3])
  p11 <- x[1] / n
  p11plus <- 1 + p11
  result <- 2 * sqrt(p11 * (1 - p11) / (n - 1)) / (p11plus * p11plus)
  if (is.null(listNames)) {
    names(result) <- NULL
  } else {
    names(result) <- paste("Sorensen-Dice disimilarity upper confidence limit ", listNames[1], ",", listNames[2], sep = "")
  }
  return(result)
}

#' @describeIn seSorensen S3 method for class "numeric"
#' @export
seSorensen.numeric <- function(x, listNames = NULL, check.table = TRUE, ...) {
  if (check.table){
    if (!nice2x2Table.numeric(x)) {
      print(x)
      stop("Inadequate table to compute the standard error")
    }
  }
  n <- sum(x[1:3])
  p11 <- x[1] / n
  p11plus <- 1 + p11
  result <- 2 * sqrt(p11 * (1 - p11) / (n - 1)) / (p11plus * p11plus)
  if (is.null(listNames)) {
    names(result) <- NULL
  } else {
    names(result) <- paste("Sorensen-Dice disimilarity upper confidence limit ", listNames[1], ",", listNames[2], sep = "")
  }
  return(result)
}

#' @describeIn seSorensen S3 method for class "character"
#' @export
seSorensen.character <- function(x, y, listNames = c("gene.list1", "gene.list2"),
                                 check.table = TRUE, ...){
  tab <- buildEnrichTable(x, y, listNames, check.table, ...)
  # Typical ... arguments:
  # geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
  # onto = onto, GOLevel = ontoLevel,
  result <- seSorensen.table(tab, listNames = listNames, check.table = FALSE)
  return(result)
}

#' @describeIn seSorensen S3 method for class "list"
#' @export
seSorensen.list <- function(x, check.table = TRUE, ...){
  numLists <- length(x)
  lstNams <- names(x)
  result <- matrix(0.0, ncol = numLists, nrow = numLists)
  for (iLst1 in 2:numLists) {
    for (iLst2 in 1:(iLst1-1)) {
      result[iLst1, iLst2] <- seSorensen.character(x[[iLst1]], x[[iLst2]],
                                                   listNames = NULL, check.table = check.table, ...)
    }
  }
  result[upper.tri(result)] <- t(result)[upper.tri(result)]
  colnames(result) <- lstNams
  rownames(result) <- lstNams
  return(result)
}
