#' Standard error of the sample Sorensen-Dice dissimilarity, asymptotic approach
#'
#' @param x either an object of class "table", "matrix" or "numeric" representing a 2x2 contingency table,
#' or a "character" (a set of gene identifiers) or "list" object. See the details section for more information.
#' @param y an object of class "character" representing a vector of gene identifiers.
#' @param n numeric. It is ignored except in the "table", "matrix" or "numeric" interfaces when argument \code{x}
#' represents relative frequencies, see the details section for more information.
#' @param check.table Boolean. If TRUE (default), argument \code{x} is checked to adequately
#' represent a 2x2 contingency table. This checking is performed by means of function
#' \code{nice2x2Table}  (only in the "table", "matrix" or "numeric" interfaces).
#' @param ... extra parameters for function \code{crossTabGOIDs4GeneLists} in package \code{equivStandardTest}.
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
#'   n_{00} \tab n_{01} \cr
#'   n_{10} \tab n_{11},
#'  }
#' }{}
#'
#'\tabular{rr}{
#' n_00 \tab n_01 \cr
#' n_10 \tab n_11,
#'}
#'
#' The subindex '00' corresponds to those GO items non enriched in both gene lists, '10' corresponds to those
#' enriched in the first list but not in the second, '01' to items non enriched in the first list
#' but enriched in the second and '11' to those GO items enriched in both lists.
#' These values must be provided in this order. In the "numeric" interface,
#' if \code{length(x) == 3} the values are interpreted as
#' \eqn{(n_{10}, n_{01}, n_{11})}{%
#' (n_10, n_01, n_11)} and otherwise as
#' \eqn{(n_{00}, n_{10}, n_{01}, n_{11})}{%
#' (n_00, n_10, n_01, n_11)}, discarding extra values if necessary.
#'
#' The result is correct, regardless the frequencies being absolute or relative but, in this second case, the value of
#' argument \code{n} must be provided and must correspond to the total number of enriched GO items, i.e., to the sum
#' of absolute frequencies
#' \eqn{n_{10} + n_{01} + n_{11}}{%
#' n_10 + n_01 + n_11.}

#'
#' @seealso \code{\link{nice2x2Table}} for checking and reformatting data,
#' \code{\link{dSorensen}} for computing the Sorensen-Dice dissimilarity,
#' \code{\link{duppSorensen}} for the upper limit of a one-sided confidence interval
#' of the dissimilarity, \code{\link{equivTestSorensen}} for an equivalence test.
#'
#' @examples
#' # Gene lists 'atlas' and 'sanger' in 'Cangenes' dataset. Table of joint enrichment
#' # of GO items in ontology BP at level 3.
#' # load(file = "atlas.sanger_BP3.Rda")
#' atlas.sanger_BP3
#' dSorensen(atlas.sanger_BP3)
#' seSorensen(atlas.sanger_BP3)
#' # To compute se for a proportions table:
#' relAtlas.sanger_BP3 <- atlas.sanger_BP3/sum(atlas.sanger_BP3)
#' seSorensen(relAtlas.sanger_BP3, n = sum(atlas.sanger_BP3[2:4]))
#'
#' library(equivStandardTest)
#' data(humanEntrezIDs)
#'
#' ?pbtGeneLists
#' # (Time consuming:)
#' seSorensen(pbtGeneLists[[2]], pbtGeneLists[[4]],
#'            listNames = names(pbtGeneLists)[c(2,4)],
#'            onto = "BP", GOLevel = 5,
#'            geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
#' # Essentially, the above code makes the same as:
#' pbtBP5.IRITD3vsKT1 <- nice2x2Table(
#'   crossTabGOIDs4GeneLists(pbtGeneLists[[2]], pbtGeneLists[[4]],
#'                           onto = "BP", GOLevel = 5,
#'                           geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db"),
#'                           listNames = names(pbtGeneLists)[c(2,4)])
#' seSorensen(pbtBP5.IRITD3vsKT1)
#'
#' # (Quite time consuming:)
#' seSorensen(pbtGeneLists,
#'            onto = "BP", GOLevel = 5,
#'            geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
#'
#' @export
seSorensen <- function(x, ...) {
  UseMethod("seSorensen")
}

#' @describeIn seSorensen S3 method for class "table"
#' @export
seSorensen.table <- function(x, n, check.table = TRUE) {
  if (check.table){
    x <- nice2x2Table(x)
  }
  if (missing(n)) {
    n <- sum(x[2:4])
    pij_samp <- x[2:4] / n
  } else {
    pij_samp <- x[2:4] / sum(x[2:4])
  }
  pNonCoincid <- pij_samp[1] + pij_samp[2]
  denom <- (2 * pij_samp[3] + pNonCoincid)^2
  t11 <- -2 * pNonCoincid / denom
  t10 <- t01 <- 2 * pij_samp[3] / denom
  sig2 <- 4 * pij_samp[3] * pNonCoincid / (denom * denom)
  return(sqrt(sig2 / n))
}

#' @describeIn seSorensen S3 method for class "matrix"
#' @export
seSorensen.matrix <- function(x, n, check.table = TRUE) {
  # x <- as.table(x)
  # seSorensen.table(x, n, check.table)
  if (check.table){
    x <- nice2x2Table(x)
  }
  if (missing(n)) {
    n <- sum(x[2:4])
    pij_samp <- x[2:4] / n
  } else {
    pij_samp <- x[2:4] / sum(x[2:4])
  }
  pNonCoincid <- pij_samp[1] + pij_samp[2]
  denom <- (2 * pij_samp[3] + pNonCoincid)^2
  t11 <- -2 * pNonCoincid / denom
  t10 <- t01 <- 2 * pij_samp[3] / denom
  sig2 <- 4 * pij_samp[3] * pNonCoincid / (denom * denom)
  return(sqrt(sig2 / n))
}

#' @describeIn seSorensen S3 method for class "numeric"
#' @export
seSorensen.numeric <- function(x, n, check.table = TRUE) {
  if (check.table){
    x <- nice2x2Table(x, n)
  }
  if (length(x) == 3) {
    if (missing(n)) {
      n <- sum(x)
      pij_samp <- x / n
    } else {
      pij_samp <- x / sum(x)
    }
  } else {
    if (missing(n)) {
      n <- sum(x[2:4])
      pij_samp <- x[2:4] / n
    } else {
      pij_samp <- x[2:4] / sum(x[2:4])
    }
  }
  pNonCoincid <- pij_samp[1] + pij_samp[2]
  denom <- (2 * pij_samp[3] + pNonCoincid)^2
  t11 <- -2 * pNonCoincid / denom
  t10 <- t01 <- 2 * pij_samp[3] / denom
  sig2 <- 4 * pij_samp[3] * pNonCoincid / (denom * denom)
  result <- sqrt(sig2 / n)
  names(result) <- NULL
  return(result)
}

#' @describeIn seSorensen S3 method for class "character"
#' @export
seSorensen.character <- function(x, y, listNames = c("gene.list1", "gene.list2"),
                                ...){
  tab <- crossTabGOIDs4GeneLists (genelist1 = x, genelist2 = y, ...)
  # Typical ... arguments:
  # geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
  # onto = onto, GOLevel = ontoLevel,
  tab <- nice2x2Table.table(tab, listNames)
  result <- seSorensen.table(tab, check.table = FALSE)
  if (is.null(listNames)) {
    names(result) <- NULL
  } else {
    names(result) <- paste("Sorensen-Dice disimilarity standard error ", listNames[1], ",", listNames[2], sep = "")
  }
  return(result)
}

#' @describeIn seSorensen S3 method for class "list"
#' @export
seSorensen.list <- function(x, ...){
  numLists <- length(x)
  lstNams <- names(x)
  result <- matrix(0.0, ncol = numLists, nrow = numLists)
  for (iLst1 in 2:numLists) {
    for (iLst2 in 1:(iLst1-1)) {
      result[iLst1, iLst2] <- seSorensen.character(x[[iLst1]], x[[iLst2]], listNames = NULL, ...)
    }
  }
  result[upper.tri(result)] <- t(result)[upper.tri(result)]
  colnames(result) <- lstNams
  rownames(result) <- lstNams
  return(result)
}
