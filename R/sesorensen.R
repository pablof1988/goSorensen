#' Standard error of the sample Sorensen-Dice dissimilarity, asymptotic approach
#'
#' @param x either an object of class "table", "matrix" or "numeric" representing a 2x2 contingency table,
#' or a "character" (a set of gene identifiers) or "list" or "tableList" object.
#' See the details section for more information.
#' @param y an object of class "character" representing a vector of gene identifiers.
#' @param check.table Boolean. If TRUE (default), argument \code{x} is checked to adequately
#' represent a 2x2 contingency table. This checking is performed by means of function
#' \code{nice2x2Table}.
#' @param ... extra parameters for function \code{buildEnrichTable}.
#'
#' @return In the "table", "matrix", "numeric" and "character" interfaces, the value of the standard error of the
#' Sorensen-Dice dissimilarity estimate.
#' In the "list" and "tableList" interfaces, the symmetric matrix of all
#' standard error dissimilarity estimates.
#'
#' @details
#' This function computes the standard error estimate of the sample Sorensen-Dice dissimilarity,
#' given a 2x2 arrangement of frequencies (either implemented as a "table", a "matrix"
#' or a "numeric" object):
#'\tabular{rr}{
#' n11 \tab n10 \cr
#' n01 \tab n00,
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
#' If \code{x} is an object of class "character", then \code{x} (and \code{y}) must represent
#' two "character" vectors of valid gene identifiers.
#' Then the standard error for the dissimilarity between lists \code{x} and \code{y} is computed,
#' after internally summarizing them as a 2x2 contingency table of joint enrichment.
#' This last operation is performed by function \code{\link{buildEnrichTable}} and "valid gene
#' identifiers" stands for the coherency of these gene identifiers with the arguments
#' \code{geneUniverse} and \code{orgPackg} of \code{buildEnrichTable}, passed by the ellipsis
#' argument \code{...} in \code{seSorensen}.
#'
#' In the "list" interface, the argument must be a list of "character" vectors, each one
#' representing a gene list (character identifiers). Then, all pairwise standard errors
#' of the dissimilarity between these gene lists are computed.
#'
#' If \code{x} is an object of class "tableList", the standard error of the Sorensen-Dice dissimilarity
#' estimate is computed over each one of these tables.
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
#' \code{\link{nice2x2Table}} for checking the validity of enrichment contingency tables,
#' \code{\link{dSorensen}} for computing the Sorensen-Dice dissimilarity,
#' \code{\link{duppSorensen}} for the upper limit of a one-sided confidence interval
#' of the dissimilarity, \code{\link{equivTestSorensen}} for an equivalence test.
#'
#' @examples
#' # Gene lists 'atlas' and 'sanger' in 'Cangenes' dataset. Table of joint enrichment
#' # of GO items in ontology BP at level 3.
#' data(tab_atlas.sanger_BP3)
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
#' #            onto = "CC", GOLevel = 5,
#' #            geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
#' # Essentially, the above code makes the same as:
#' # tab.IRITD3vsKT1 <- buildEnrichTable(pbtGeneLists[[2]], pbtGeneLists[[4]],
#' #                                     onto = "CC", GOLevel = 5,
#' #                                     geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
#' # tab.IRITD3vsKT1
#' # seSorensen(tab.IRITD3vsKT1)
#'
#' # All pairwise standard errors (quite time consuming):
#' # seSorensen(pbtGeneLists,
#' #            onto = "CC", GOLevel = 5,
#' #            geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
#' @md

#' @export
seSorensen <- function(x, ...) {
  UseMethod("seSorensen")
}

#' @describeIn seSorensen S3 method for class "table"
#' @export
seSorensen.table <- function(x, check.table = TRUE, ...) {
  stopifnot("Argument 'check.table' must be logical" = is.logical(check.table))
  if (check.table){
    nice2x2Table.table(x)
  }
  n <- sum(x[seq.int(1, 3)])
  p11 <- x[1] / n
  p11plus <- 1 + p11
  result <- 2 * sqrt(p11 * (1 - p11) / (n - 1)) / (p11plus * p11plus)
  return(result)
}

#' @describeIn seSorensen S3 method for class "matrix"
#' @export
seSorensen.matrix <- function(x, check.table = TRUE, ...) {
  stopifnot("Argument 'check.table' must be logical" = is.logical(check.table))
  if (check.table){
    nice2x2Table.matrix(x)
  }
  n <- sum(x[seq.int(1, 3)])
  p11 <- x[1] / n
  p11plus <- 1 + p11
  result <- 2 * sqrt(p11 * (1 - p11) / (n - 1)) / (p11plus * p11plus)
  return(result)
}

#' @describeIn seSorensen S3 method for class "numeric"
#' @export
seSorensen.numeric <- function(x, check.table = TRUE, ...) {
  stopifnot("Argument 'check.table' must be logical" = is.logical(check.table))
  if (check.table){
    nice2x2Table.numeric(x)
  }
  n <- sum(x[seq.int(1, 3)])
  p11 <- x[1] / n
  p11plus <- 1 + p11
  result <- 2 * sqrt(p11 * (1 - p11) / (n - 1)) / (p11plus * p11plus)
  return(result)
}

#' @describeIn seSorensen S3 method for class "character"
#' @export
seSorensen.character <- function(x, y, check.table = TRUE, ...){
  tab <- buildEnrichTable(x, y, check.table = check.table, ...)
  # Typical ... arguments:
  # geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
  # onto = onto, GOLevel = ontoLevel,
  result <- seSorensen.table(tab, check.table = FALSE)
  return(result)
}

#' @describeIn seSorensen S3 method for class "list"
#' @export
seSorensen.list <- function(x, check.table = TRUE, ...){
  tabs <- buildEnrichTable(x, check.table = check.table, ...)
  return(seSorensen.tableList(tabs, check.table = FALSE))
}

#' @describeIn seSorensen S3 method for class "tableList"
#' @export
seSorensen.tableList <- function(x, check.table = TRUE, ...){
  numLists <- length(x) + 1
  lstNams <- c(names(x[[1]]), names(x))
  uTri <- sapply(seq.int(1, numLists - 1), function(iLst1) {
    aux <- vapply(seq.int(1, iLst1), function(iLst2) {
      seSorensen.table(x[[iLst1]][[iLst2]],
                       check.table = check.table)
    }, FUN.VALUE = 0.0)
    return(aux)
  })
  result <- matrix(0.0, ncol = numLists, nrow = numLists)
  result[upper.tri(result)] <- unlist(uTri)
  result[lower.tri(result)] <- t(result)[lower.tri(result)]
  colnames(result) <- lstNams
  rownames(result) <- lstNams
  return(result)
}
