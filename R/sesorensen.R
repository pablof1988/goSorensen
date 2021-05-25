#' Standard error of the sample Sorensen-Dice dissimilarity, asymptotic approach
#'
#' @param x either an object of class "table", "matrix" or "numeric" representing a 2x2 contingency table,
#' or a "character" (a set of gene identifiers) or "list" object. See the details section for more information.
#' @param y an object of class "character" representing a vector of gene identifiers.
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
#' It is based on an asymptotic approach.
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
#' seSorensen(atlas.sanger_BP3/sum(atlas.sanger_BP3)) / sqrt(sum(atlas.sanger_BP3))
#' # In general:
#' # seSorensen(pij) / sqrt(n)
#'
#' library(equivStandardTest)
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
seSorensen.table <- function(x, n = sum(x), check.table = TRUE) {
  if (check.table){
    x <- nice2x2Table(x)
  }
  pij_samp <- x[2:4] / n
  denom <- (2 * pij_samp[3] + pij_samp[1] + pij_samp[2])^2
  t11 <- -2 * (pij_samp[1] + pij_samp[2]) / denom
  t21 <- t12 <- 2 * pij_samp[3] / denom
  sig_cuad <- (t11^2 * pij_samp[3] + t12^2 * pij_samp[1] + t21^2 * pij_samp[2]) -
    (t11 * pij_samp[3] + t12 * pij_samp[1] + t21 * pij_samp[2])^2
  return(sqrt(sig_cuad / n))
}

#' @describeIn seSorensen S3 method for class "matrix"
#' @export
seSorensen.matrix <- function(x, n = sum(x), check.table = TRUE) {
  x <- as.table(x)
  seSorensen.table(x, n, check.table)
}

#' @describeIn seSorensen S3 method for class "numeric"
#' @export
seSorensen.numeric <- function(x, n, check.table = TRUE) {
  if (check.table){
    x <- nice2x2Table(x, n)
  }
  n <- sum(x)
  pij_samp <- x[2:4] / n
  denom <- (2 * pij_samp[3] + pij_samp[1] + pij_samp[2])^2
  t11 <- -2 * (pij_samp[1] + pij_samp[2]) / denom
  t21 <- t12 <- 2 * pij_samp[3] / denom
  sig_cuad <- (t11^2 * pij_samp[3] + t12^2 * pij_samp[1] + t21^2 * pij_samp[2]) -
    (t11 * pij_samp[3] + t12 * pij_samp[1] + t21 * pij_samp[2])^2
  result <- sqrt(sig_cuad / n)
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
