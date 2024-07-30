#' Computation of the Sorensen-Dice dissimilarity
#' @param x either an object of class "table", "matrix" or "numeric" representing a 2x2 contingency table,
#' or a "character" vector (a set of gene identifiers) or "list" or "tableList" object.
#' See the details section for more information.
#' @param y an object of class "character" representing a vector of valid gene identifiers (e.g., ENTREZ).
#' @param check.table Boolean. If TRUE (default), argument \code{x} is checked to adequately
#' represent a 2x2 contingency table, by means of function \code{nice2x2Table}.
#' @param ... extra parameters for function \code{buildEnrichTable}.
#'
#' @return In the "table", "matrix", "numeric" and "character" interfaces, the value of the
#' Sorensen-Dice dissimilarity. In the "list" and "tableList" interfaces, the symmetric matrix
#' of all pairwise Sorensen-Dice dissimilarities.
#'
#' @details
#' Given a 2x2 arrangement of frequencies (either implemented as a "table", a "matrix"
#' or a "numeric" object):
#'
#'\tabular{rr}{
#' \eqn{n_{11}} \tab \eqn{n_{10}} \cr
#' \eqn{n_{01}} \tab \eqn{n_{00}},
#'}
#'
#' this function computes the Sorensen-Dice dissimilarity
#' \deqn{\frac{n_{10} + n_{01}}{2 n_{11} + n_{10} + n_{01}}.}{{%
#' n_10 + n_01}/{2 n_11 + n_10 + n_01}.}
#'
#' The subindex '11' corresponds to those
#' GO terms enriched in both lists, '01' to terms enriched in the second list but not in the first one,
#' '10' to terms enriched in the first list but not enriched in the second one and '00' corresponds
#' to those GO terms non enriched in both gene lists, i.e., to the double negatives, a value which
#' is ignored in the computations.
#'
#' In the "numeric" interface, if \code{length(x) >= 3}, the values are interpreted
#' as
#' \eqn{(n_{11}, n_{01}, n_{10}, n_{00})}{%
#' (n_11, n_01, n_10, n_00)}, always in this order and discarding extra values if necessary.
#' The result is correct, regardless the frequencies being absolute or relative.
#'
#' If \code{x} is an object of class "character", then \code{x} (and \code{y}) must represent
#' two "character" vectors of valid gene identifiers (e.g., ENTREZ).
#' Then the dissimilarity between lists \code{x} and \code{y} is computed,
#' after internally summarizing them as a 2x2 contingency table of joint enrichment.
#' This last operation is performed by function \code{\link{buildEnrichTable}} and "valid gene
#' identifiers (e.g., ENTREZ)" stands for the coherency of these gene identifiers with the arguments
#' \code{geneUniverse} and \code{orgPackg} of \code{buildEnrichTable}, passed by the ellipsis
#' argument \code{...} in \code{dSorensen}.
#'
#' If \code{x} is an object of class "list", the argument must be a list of "character" vectors, each one
#' representing a gene list (character identifiers). Then, all pairwise dissimilarities between these
#' gene lists are computed.
#'
#' If \code{x} is an object of class "tableList", the Sorensen-Dice dissimilarity is computed
#' over each one of these tables.
#' Given k gene lists (i.e. "character" vectors of gene identifiers) l1, l2, ..., lk,
#' an object of class "tableList" (typically constructed by a call to function
#' \code{\link{buildEnrichTable}}) is a list of lists of
#' contingency tables t(i,j) generated from each pair of gene lists i and j, with the
#' following structure:
#'
#' $l2
#'
#' $l2$l1$t(2,1)
#'
#' $l3
#'
#' $l3$l1$t(3,1), $l3$l2$t(3,2)
#'
#' ...
#'
#' $lk
#'
#' $lk$l1$t(k,1), $lk$l2$t(k,2), ..., $lk$l(k-1)t(k,k-1)
#'
#' @seealso
#' \code{\link{buildEnrichTable}} for constructing contingency tables of mutual
#' enrichment,
#' \code{\link{nice2x2Table}} for checking contingency tables validity,
#' \code{\link{seSorensen}} for computing the standard error of the dissimilarity,
#' \code{\link{duppSorensen}} for the upper limit of a one-sided confidence interval
#' of the dissimilarity, \code{\link{equivTestSorensen}} for an equivalence test.
#'
#' @examples
#' # Gene lists 'atlas' and 'sanger' in 'allOncoGeneLists' dataset. Table of joint enrichment
#' # of GO terms in ontology BP at level 3.
#' data(tab_atlas.sanger_BP3)
#' tab_atlas.sanger_BP3
#' ?tab_atlas.sanger_BP3
#' dSorensen(tab_atlas.sanger_BP3)
#'
#' # Table represented as a vector:
#' conti4 <- c(56, 1, 30, 471)
#' dSorensen(conti4)
#' # or as a plain matrix:
#' dSorensen(matrix(conti4, nrow = 2))
#'
#' # This function is also appropriate for proportions:
#' dSorensen(conti4 / sum(conti4))
#'
#' conti3 <- c(56, 1, 30)
#' dSorensen(conti3)
#'
#' # Sorensen-Dice dissimilarity from scratch, directly from two gene lists:
#' # (These examples may be considerably time consuming due to many enrichment
#' # tests to build the contingency tables of joint enrichment)
#' # data(allOncoGeneLists)
#' # ?allOncoGeneLists
#' 
#' # Obtaining ENTREZ identifiers for the gene universe of humans:
#' # library(org.Hs.eg.db)
#' # humanEntrezIDs <- keys(org.Hs.eg.db, keytype = "ENTREZID")
#' 
#' # (Time consuming, building the table requires many enrichment tests:)
#' # dSorensen(allOncoGeneLists$atlas, allOncoGeneLists$sanger,
#' #           onto = "BP", GOLevel = 3,
#' #           geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
#' 
#' # Essentially, the above code makes the same as:
#' # tab_atlas.sanger_BP3 <- buildEnrichTable(allOncoGeneLists$atlas, allOncoGeneLists$sanger,
#' #                                     onto = "BP", GOLevel = 3,
#' #                                     geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
#' # dSorensen(tab_atlas.sanger_BP3)
#' # (Quite time consuming, all pairwise dissimilarities:)
#' # dSorensen(allOncoGeneLists,
#' #           onto = "BP", GOLevel = 3,
#' #           geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
#' @md

#' @export
dSorensen <- function(x, ...) {
  UseMethod("dSorensen")
}

#' @describeIn dSorensen S3 method for class "table"
#' @export
dSorensen.table <- function(x, check.table = TRUE, ...){
  stopifnot("Argument 'check.table' must be logical" = is.logical(check.table))
  if (check.table){
    nice2x2Table.table(x)
  }
  result <- (x[2] + x[3]) / (2 * x[1] + x[2] + x[3])
  return(result)
}

#' @describeIn dSorensen S3 method for class "matrix"
#' @export
dSorensen.matrix <- function(x, check.table = TRUE, ...){
  stopifnot("Argument 'check.table' must be logical" = is.logical(check.table))
  if (check.table){
    nice2x2Table.matrix(x)
  }
  result <- (x[2] + x[3]) / (2 * x[1] + x[2] + x[3])
  return(result)
}

#' @describeIn dSorensen S3 method for class "numeric"
#' @export
dSorensen.numeric <- function(x, check.table = TRUE, ...){
  stopifnot("Argument 'check.table' must be logical" = is.logical(check.table))
  if (check.table){
    nice2x2Table.numeric(x)
  }
  result <- (x[2] + x[3]) / (2 * x[1] + x[2] + x[3])
  return(result)
}

#' @describeIn dSorensen S3 method for class "character"
#' @export
dSorensen.character <- function(x, y, check.table = TRUE, ...){
  tab <- buildEnrichTable(x, y, check.table = check.table, ...)
  # Typical ... arguments:
  # geneUniverse=humanEntrezIDs, orgPackg="org.Hs.eg.db",
  # onto = onto, GOLevel = ontoLevel,
  result <- (tab[2] + tab[3]) / (2 * tab[1] + tab[2] + tab[3])
  return(result)
}

#' @describeIn dSorensen S3 method for class "list"
#' @export
dSorensen.list <- function(x, check.table = TRUE, ...)
{
  tabs <- buildEnrichTable(x, check.table = check.table, ...)
  return(dSorensen.tableList(tabs, check.table = FALSE))
}

#' @describeIn dSorensen S3 method for class "tableList"
#' @export
dSorensen.tableList <- function(x, check.table = TRUE, ...){
  numLists <- length(x) + 1
  lstNams <- c(names(x[[1]]), names(x))
  uTri <- sapply(seq.int(1, numLists - 1), function(iLst1) {
    vapply(seq.int(1, iLst1), function(iLst2) {
      dSorensen.table(x[[iLst1]][[iLst2]],
                      check.table = check.table)
    }, FUN.VALUE = 0.0)
  })
  result <- matrix(0.0, ncol = numLists, nrow = numLists)
  result[upper.tri(result)] <- unlist(uTri)
  result[lower.tri(result)] <- t(result)[lower.tri(result)]
  colnames(result) <- lstNams
  rownames(result) <- lstNams
  return(result)
}
