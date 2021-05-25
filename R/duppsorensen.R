#' Upper limit of a one-sided confidence interval (0, dUpp] for the Sorensen-Dice dissimilarity
#'
#' @param x either an object of class "table", "matrix" or "numeric" representing a 2x2 contingency table,
#' or a "character" (a set of gene identifiers) or "list" object. See the details section for more information.
#' @param y an object of class "character" representing a vector of gene identifiers.
#' @param n total number of cases (GO nodes) in the contingency table. Only required
#' (sometimes) on the "numeric" interface, see the details section.
#' @param dis Sorensen-Dice dissimilarity value. Only required to speed computations if this value
#' is known in advance.
#' @param se standard error estimate of the sample dissimilarity. Only required to speed computations
#' if this value is known in advance.
#' @param conf.level confidence level of the one-sided confidence interval, a value between 0 and 1.
#' @param z.conf.level standard normal distribution quantile at the \code{conf.level} value. Only
#' required to speed computations if this value is known in advance. Then, the parameter \code{conf.level}
#' is ignored.
#' @param check.table Boolean. If TRUE (default), argument \code{x} is checked to adequately
#' represent a 2x2 contingency table. This checking is performed by means of function
#' \code{nice2x2Table} (only in the "table", "matrix" or "numeric" interfaces).
#'
#' @return In the "table", "matrix", "numeric" and "character" interfaces, the value of the Upper limit of the confidence
#' interval for the Sorensen-Dice dissimilarity. In the "list" interface, the symmetric matrix of all pairwise upper limits.
#'
#' @details
#' Arguments \code{dis}, \code{se} and \code{z.conf.level} are not required. If known in advance (e.g., as
#' a consequence of previous computations with the same data), providing its value may speed the computations.
#' If \code{x} is an object of class "character", it must represent a list of gene identifiers. Then the
#' confidence interval for the dissimilarity between lists \code{x} and \code{y} is computed, after summarizing
#' these gene lists as a 2x2 contingency table of joint enrichment.
#' In the "list" interface, the argument must be a list of "character" vectors, each one
#' representing a gene list (character identifiers). Then, all pairwise upper limits of the dissimilarity between
#' these gene lists are computed.
#'
#' @seealso \code{\link{nice2x2Table}} for checking and reformatting data,
#' \code{\link{dSorensen}} for computing the Sorensen-Dice dissimilarity,
#' \code{\link{seSorensen}} for computing the standard error of the dissimilarity,
#' \code{\link{equivTestSorensen}} for an equivalence test.
#'
#' @examples
#' # Gene lists 'atlas' and 'sanger' in 'Cangenes' dataset. Table of joint enrichment
#' # of GO items in ontology BP at level 3.
#' tab_atlas.sanger_BP3
#' dSorensen(tab_atlas.sanger_BP3)
#' seSorensen(tab_atlas.sanger_BP3)
#' duppSorensen(tab_atlas.sanger_BP3)
#'
#' library(equivStandardTest)
#' ?pbtGeneLists
#' # (Time consuming:)
#' duppSorensen(pbtGeneLists[[2]], pbtGeneLists[[4]],
#'              listNames = names(pbtGeneLists)[c(2,4)],
#'              onto = "BP", GOLevel = 5,
#'              geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
#' # Essentially, the above code makes the same as:
#' pbtBP5.IRITD3vsKT1 <- nice2x2Table(
#'   crossTabGOIDs4GeneLists(pbtGeneLists[[2]], pbtGeneLists[[4]],
#'                           onto = "BP", GOLevel = 5,
#'                           geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db"),
#'                           listNames = names(pbtGeneLists)[c(2,4)])
#' duppSorensen(pbtBP5.IRITD3vsKT1)
#' # (Quite time consuming:)
#' duppSorensen(pbtGeneLists,
#'              onto = "BP", GOLevel = 5,
#'              geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
#'
#' @export
duppSorensen <- function(x, ...) {
  UseMethod("duppSorensen")
}

#' @describeIn duppSorensen S3 method for class "table"
#' @export
duppSorensen.table <- function(x, n = sum(x), dis = dSorensen.table(x),
                               se = seSorensen.table(x, n, FALSE),
                               conf.level = 0.95, z.conf.level = qnorm(conf.level),
                               check.table = TRUE){
  if (check.table){
    x <- nice2x2Table(x)
  }
  return(dis + z.conf.level * se)
}

#' @describeIn duppSorensen S3 method for class "matrix"
#' @export
duppSorensen.matrix <- function(x, n = sum(x), dis = dSorensen.matrix(x),
                                se = seSorensen.matrix(x, n, FALSE),
                                conf.level = 0.95, z.conf.level = qnorm(conf.level),
                                check.table = TRUE){
  if (check.table){
    x <- nice2x2Table(x)
  }
  return(dis + z.conf.level * se)
}

#' @describeIn duppSorensen S3 method for class "numeric"
#' @export
duppSorensen.numeric <- function(x, n, dis = dSorensen.numeric(x, n, FALSE),
                                 se = seSorensen.numeric(x, n, FALSE),
                                 conf.level = 0.95, z.conf.level = qnorm(conf.level),
                                 check.table = TRUE){
  if (check.table){
    x <- nice2x2Table(x, n)
  }
  return(dis + z.conf.level * se)
}

#' @describeIn duppSorensen S3 method for class "character"
#' @export
duppSorensen.character <- function(x, y, conf.level = 0.95,
                                   dis, se, z.conf.level = qnorm(conf.level),
                                   listNames = c("gene.list1", "gene.list2"), check.table = TRUE,
                                   ...){
  tab <- crossTabGOIDs4GeneLists (genelist1 = x, genelist2 = y, ...)
  # Typical ... arguments:
  # geneUniverse=humanEntrezIDs, orgPackg="org.Hs.eg.db",
  # onto = onto, GOLevel = ontoLevel,
  tab <- nice2x2Table.table(tab, listNames)
  result <- duppSorensen.table(tab,
                               dis = if (missing(dis)) dSorensen.table(tab, check.table = FALSE) else dis,
                               se = if (missing(se)) seSorensen.table(tab, check.table = FALSE) else se,
                               z.conf.level = if (missing(z.conf.level)) qnorm(conf.level) else z.conf.level,
                               conf.level = 0.95)
  if (is.null(listNames)) {
    names(result) <- NULL
  } else {
    names(result) <- paste("Sorensen-Dice disimilarity upper confidence limit ", listNames[1], ",", listNames[2], sep = "")
  }
  return(result)
}

#' @describeIn duppSorensen S3 method for class "list"
#' @export
duppSorensen.list <- function(x, ...){
  numLists <- length(x)
  lstNams <- names(x)
  result <- matrix(0.0, ncol = numLists, nrow = numLists)
  for (iLst1 in 2:numLists) {
    for (iLst2 in 1:(iLst1-1)) {
      result[iLst1, iLst2] <- duppSorensen.character(x[[iLst1]], x[[iLst2]], listNames = NULL, ...)
    }
  }
  colnames(result) <- lstNams
  rownames(result) <- lstNams
  return(result)
}
