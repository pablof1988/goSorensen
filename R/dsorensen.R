#' Computation of the Sorensen-Dice dissimilarity
#' @param x either an object of class "table", "matrix" or "numeric" representing a 2x2 contingency table,
#' or a "character" (a set of gene identifiers) or "list" object. See the details section for more information.
#' @param y an object of class "character" representing a vector of gene identifiers.
#' @param check.table Boolean. If TRUE (default), argument \code{x} is checked to adequately
#' represent a 2x2 contingency table, by means of function \code{nice2x2Table}  (only in the "table", "matrix"
#' or "numeric" interfaces).
#' @param ... extra parameters for function \code{crossTabGOIDs4GeneLists} in package \code{equivStandardTest}.
#'
#' @return In the "table", "matrix", "numeric" and "character" interfaces, the value of the Sorensen-Dice
#' dissimilarity. In the "list" interface, the symmetric matrix of all pairwise dissimilarities.
#'
#' @details
#' Given a 2x2 arrangement of frequencies (either implemented as a "table", a "matrix"
#' or a "numeric" object):
#' \tabular{rr}{
#'  n_{00} \tab n_{01} \cr
#'  n_{10} \tab n_{11}
#' }
#' this function computes the Sorensen-Dice dissimilarity
#' \eqn{(n_{10} + n_{01}) / (2 n_{11} + n_{10} + n_{01})}.
#' The result is correct, regardless the frequencies are absolute or relative.
#' If \code{x} is an object of class "character", it must represent a list of gene identifiers. Then the
#' dissimilarity between lists \code{x} and \code{y} is computed, after summarizing these gene lists
#' as a 2x2 contingency table of joint enrichment.
#' In the "list" interface, the argument must be a list of "character" vectors, each one
#' representing a gene list (character identifiers). Then, all pairwise dissimilarities between these
#' gene lists are computed.
#'
#' @seealso \code{\link{nice2x2Table}} for checking and reformatting data,
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
#' # Badly formed table:
#' conti <- as.table(matrix(c(501, 27, 36, 12, 43, 15, 0, 0, 0),
#'   nrow = 3, ncol = 3,
#'   dimnames = list(c("a1","a2","a3"),
#'                   c("b1", "b2","b3"))))
#' dSorensen(conti)
#' dSorensen(conti, check.table = FALSE)
#' # Wrong value!
#' # Correct:
#' dSorensen(nice2x2Table(conti), check.table = FALSE)
#' conti2 <- conti[1,1:min(2,ncol(conti)), drop = FALSE]
#' dSorensen(conti2)
#' dSorensen(matrix(c(210, 12), ncol = 2, nrow = 1))
#' conti4 <- c(1439, 32, 21, 81)
#' dSorensen(conti4)
#' # This function is also appropriate for proportions:
#' dSorensen(conti4 / sum(conti4))
#' dSorensen(matrix(conti4, nrow = 2))
#' conti5 <- c(32, 21, 81)
#' try(dSorensen(conti5))
#' dSorensen(conti5, n = 1573)
#' dSorensen(conti5 / 1573, n = 1)
#' dSorensen(conti5, n = 1000)
#' dSorensen(conti5 / 1000, n = 1)
#' try(dSorensen(conti5, n = 10))
#'
#' library(equivStandardTest)
#' ?pbtGeneLists
#' # (Time consuming:)
#' dSorensen(pbtGeneLists[[2]], pbtGeneLists[[4]],
#'           listNames = names(pbtGeneLists)[c(2,4)],
#'           onto = "BP", GOLevel = 5,
#'           geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
#' # Essentially, the above code makes the same as:
#' pbtBP5.IRITD3vsKT1 <- nice2x2Table(
#' crossTabGOIDs4GeneLists(pbtGeneLists[[2]], pbtGeneLists[[4]],
#'                         onto = "BP", GOLevel = 5,
#'                         geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db"),
#' listNames = names(pbtGeneLists)[c(2,4)])
#' dSorensen(pbtBP5.IRITD3vsKT1)
#' # (Quite time consuming:)
#' dSorensen(pbtGeneLists,
#'   onto = "BP", GOLevel = 5,
#'   geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")



#' @export
dSorensen <- function(x, ...) {
  UseMethod("dSorensen")
}

#' @describeIn dSorensen S3 method for class "table"
#' @export
dSorensen.table <- function(x, check.table = TRUE){
  if (check.table){
    x <- nice2x2Table(x)
  }
  return((x[2] + x[3]) / (2 * x[4] + x[2] + x[3]))
}

#' @describeIn dSorensen S3 method for class "matrix"
#' @export
dSorensen.matrix <- function(x, check.table = TRUE){
  if (check.table){
    x <- nice2x2Table(x)
  }
  return((x[2] + x[3]) / (2 * x[4] + x[2] + x[3]))
}

#' @describeIn dSorensen S3 method for class "numeric"
#' @export
dSorensen.numeric <- function(x, n, check.table = TRUE){
  if (check.table){
    x <- nice2x2Table(x, n)
  }
  result <- (x[2] + x[3]) / (2 * x[4] + x[2] + x[3])
  names(result) <- NULL
  return(result)
}

#' @describeIn dSorensen S3 method for class "character"
#' @export
dSorensen.character <- function(x, y, listNames = c("gene.list1", "gene.list2"),
                                ...){
  tab <- crossTabGOIDs4GeneLists (genelist1 = x, genelist2 = y, ...)
  # Typical ... arguments:
  # geneUniverse=humanEntrezIDs, orgPackg="org.Hs.eg.db",
  # onto = onto, GOLevel = ontoLevel,
  tab <- nice2x2Table.table(tab, listNames)
  result <- (tab[2] + tab[3]) / (2 * tab[4] + tab[2] + tab[3])
  if (is.null(listNames)) {
    names(result) <- NULL
  } else {
    names(result) <- paste("Sorensen-Dice disimilarity ", listNames[1], ",", listNames[2], sep = "")
  }
  return(result)
}

#' @describeIn dSorensen S3 method for class "list"
#' @export
dSorensen.list <- function(x, ...){
  numLists <- length(x)
  lstNams <- names(x)
  result <- matrix(0.0, ncol = numLists, nrow = numLists)
  for (iLst1 in 2:numLists) {
    for (iLst2 in 1:(iLst1-1)) {
      result[iLst1, iLst2] <- dSorensen.character(x[[iLst1]], x[[iLst2]], listNames = NULL, ...)
    }
  }
  result[upper.tri(result)] <- t(result)[upper.tri(result)]
  colnames(result) <- lstNams
  rownames(result) <- lstNams
  return(result)
}
