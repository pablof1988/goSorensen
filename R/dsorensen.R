#' Computation of the Sorensen-Dice dissimilarity
#' @param x either an object of class "table", "matrix" or "numeric" representing a 2x2 contingency table,
#' or a "character" (a set of gene identifiers) or "list" object. See the details section for more information.
#' @param y an object of class "character" representing a vector of gene identifiers.
# @param n total number of enriched GO items. Only required (sometimes) on the "numeric" interface,
# see the details section.
#' @param check.table Boolean. If TRUE (default), argument \code{x} is checked to adequately
#' represent a 2x2 contingency table, by means of function \code{nice2x2Table}.
#' @param ... extra parameters for function \code{crossTabGOIDs4GeneLists} in package \code{equivStandardTest}.
#'
#' @return In the "table", "matrix", "numeric" and "character" interfaces, the value of the Sorensen-Dice
#' dissimilarity. In the "list" interface, the symmetric matrix of all pairwise dissimilarities.
#'
#' @details
#' Given a 2x2 arrangement of frequencies (either implemented as a "table", a "matrix"
#' or a "numeric" object):
#' \deqn{
#'  \tabular{rr}{
#'   n_{11} \tab n_{01} \cr
#'   n_{10} \tab n_{00},
#'  }
#' }{}
#'
#'\tabular{rr}{
#' n_11 \tab n_01 \cr
#' n_10 \tab n_00,
#'}
#'
#' this function computes the Sorensen-Dice dissimilarity
#' \deqn{\frac{n_{10} + n_{01}}{2 n_{11} + n_{10} + n_{01}}.}{{%
#' n_10 + n_01}/{2 n_11 + n_10 + n_01}.}
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
#' (n_11, n_01, n_10)}, always in this order and discarding extra values if necessary.
#' The result is correct, regardless the frequencies are absolute or relative.
# If \code{1 <= length(x) <= 2}, \code{x[1]} must correspond to \eqn{n_{11}}{%
# n_11}, and a possible extra second value is ignored.
# Then, if \code{X[1]} is an absolute frequency, the value of argument \code{n} must
# be provided and must correspond to
# \eqn{(n_{11} + n_{01} + n_{10})}{%
# (n_11 + n_01 + n_10)}. To provide \code{n} is not necessary if \code{X[1]} corresponds to a
# relative frequency
#'
#' If \code{x} is an object of class "character", it must represent a list of gene identifiers. Then the
#' dissimilarity between lists \code{x} and \code{y} is computed, after summarizing these gene lists
#' as a 2x2 contingency table of joint enrichment.
#'
#' In the "list" interface, the argument must be a list of "character" vectors, each one
#' representing a gene list (character identifiers). Then, all pairwise dissimilarities between these
#' gene lists are computed.
#'
#' @seealso \code{\link{nice2x2Table}} for checking data,
#' \code{\link{seSorensen}} for computing the standard error of the dissimilarity,
#' \code{\link{duppSorensen}} for the upper limit of a one-sided confidence interval
#' of the dissimilarity, \code{\link{equivTestSorensen}} for an equivalence test.
#'
#' @examples
#' # Gene lists 'atlas' and 'sanger' in 'Cangenes' dataset. Table of joint enrichment
#' # of GO items in ontology BP at level 3.
#' ?tab_atlas.sanger_BP3
#' dSorensen(tab_atlas.sanger_BP3)
#'
#' try{dSorensen(matrix(c(210, 12), ncol = 2, nrow = 1))}
#' conti4 <- c(32, 21, 81, 1439)
#' dSorensen(conti4)
#' # This function is also appropriate for proportions:
#' dSorensen(conti4 / sum(conti4))
#' dSorensen(matrix(conti4, nrow = 2))
#' conti5 <- c(32, 21, 81)
#' dSorensen(conti5)
#'
#' library(equivStandardTest)
#' data(humanEntrezIDs)
#'
#' ?pbtGeneLists
#' # (Time consuming:)
#' dSorensen(pbtGeneLists[[2]], pbtGeneLists[[4]],
#'           listNames = names(pbtGeneLists)[c(2,4)],
#'           onto = "BP", GOLevel = 5,
#'           geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
#' # Essentially, the above code makes the same as:
#' pbtBP5.IRITD3vsKT1 <- buildEnrichTable(pbtGeneLists[[2]], pbtGeneLists[[4]],
#'                           onto = "BP", GOLevel = 5,
#'                           geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
#' dSorensen(pbtBP5.IRITD3vsKT1)
#' # (Quite time consuming:)
#' dSorensen(pbtGeneLists,
#'           onto = "BP", GOLevel = 5,
#'           geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")

#' @export
dSorensen <- function(x, ...) {
  UseMethod("dSorensen")
}

#' @describeIn dSorensen S3 method for class "table"
#' @export
dSorensen.table <- function(x, check.table = TRUE){
  if (check.table){
    if (!nice2x2Table.table(x)) {
      print(x)
      stop("Inadequate contingency table")}
  }
  return((x[2] + x[3]) / (2 * x[1] + x[2] + x[3]))
}

#' @describeIn dSorensen S3 method for class "matrix"
#' @export
dSorensen.matrix <- function(x, check.table = TRUE){
  if (check.table){
    if (!nice2x2Table.matrix(x)) {
      print(x)
      stop("Inadequate contingency table")}
  }
  return((x[2] + x[3]) / (2 * x[1] + x[2] + x[3]))
}

#' @describeIn dSorensen S3 method for class "numeric"
#' @export
dSorensen.numeric <- function(x, n, check.table = TRUE){
  if (check.table){
    if (!nice2x2Table.numeric(x, n)) {
      print(x)
      stop("Inadequate contingency table")}
  }
  # if (length(x) >= 3) {
  #   result <- (x[2] + x[3]) / (2 * x[1] + x[2] + x[3])
  # } else {
  #   if (x[1] < 1) {
  #     result <- (1 - x[1]) / (1 + x[1])
  #   } else {
  #     if (missing(n)) {
  #       stop("Argument 'n' required for numeric 'x' of length <= 2 representing absolute frequencies")
  #     }
  #     if (x[1] > n) {
  #       stop("Invalid value for argument 'n', it must be n_11 + n_10 + n_01")
  #     }
  #     p11 <- x[1] / n
  #     result <- (1 - p11) / (1 + p11)
  #   }
  # }
  result <- result <- (x[2] + x[3]) / (2 * x[1] + x[2] + x[3])
  names(result) <- NULL
  return(result)
}

#' @describeIn dSorensen S3 method for class "character"
#' @export
dSorensen.character <- function(x, y,
                                listNames = c("gene.list1", "gene.list2"),
                                check.table = TRUE,
                                ...){
  tab <- buildEnrichTable(x, y, listNames, check.table, ...)
  # Typical ... arguments:
  # geneUniverse=humanEntrezIDs, orgPackg="org.Hs.eg.db",
  # onto = onto, GOLevel = ontoLevel,
  # tab <- nice2x2Table.table(tab, listNames)
  result <- (tab[2] + tab[3]) / (2 * tab[1] + tab[2] + tab[3])
  if (is.null(listNames)) {
    names(result) <- NULL
  } else {
    names(result) <- paste("Sorensen-Dice disimilarity upper confidence limit ", listNames[1], ",", listNames[2], sep = "")
  }
  return(result)
}

#' @describeIn dSorensen S3 method for class "list"
#' @export
dSorensen.list <- function(x, check.table = TRUE, ...){
  numLists <- length(x)
  lstNams <- names(x)
  result <- matrix(0.0, ncol = numLists, nrow = numLists)
  for (iLst1 in 2:numLists) {
    for (iLst2 in 1:(iLst1-1)) {
      result[iLst1, iLst2] <- dSorensen.character(x[[iLst1]], x[[iLst2]], check.table = check.table,
                                                  listNames = NULL, ...)
    }
  }
  result[upper.tri(result)] <- t(result)[upper.tri(result)]
  colnames(result) <- lstNams
  rownames(result) <- lstNams
  return(result)
}
