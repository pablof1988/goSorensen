#' Upper limit of a one-sided confidence interval (0, dUpp] for the Sorensen-Dice dissimilarity
#'
#' @param x either an object of class "table", "matrix" or "numeric" representing a 2x2 contingency table,
#' or a "character" (a set of gene identifiers) or "list" object. See the details section for more information.
#' @param y an object of class "character" representing a vector of gene identifiers.
#' @param dis Sorensen-Dice dissimilarity value. Only required to speed computations if this value
#' is known in advance.
#' @param se standard error estimate of the sample dissimilarity. Only required to speed computations
#' if this value is known in advance.
#' @param conf.level confidence level of the one-sided confidence interval, a numeric value between 0 and 1.
#' @param z.conf.level standard normal (or bootstrap, see arguments below) distribution quantile at the
#' \code{1 - conf.level} value. Only required to speed computations if this value is known in advance.
#' Then, the argument \code{conf.level} is ignored.
#' @param boot boolean. If TRUE, \code{z.conf.level} is computed by means of a bootstrap
#' approach instead of the asymptotic normal approach. Defaults to FALSE.
#' @param nboot numeric, number of bootstrap replicates. Ignored if \code{boot == FALSE}. Defaults to 10000.
#' @param check.table Boolean. If TRUE (default), argument \code{x} is checked to adequately
#' represent a 2x2 contingency table. This checking is performed by means of function
#' \code{nice2x2Table}.
#' @param listNames character(2), names of both gene lists to be compared (only for the character interface)
#' @param ... additional arguments for function \code{buildEnrichTable}.
#'
#' @return In the "table", "matrix", "numeric" and "character" interfaces, the value of the Upper limit of the confidence
#' interval for the Sorensen-Dice dissimilarity. In the "list" interface, the symmetric matrix of all pairwise upper limits.
#'
#' @details
#' This function computes the upper limit of a one-sided confidence interval for the Sorensen-Dice dissimilarity,
#' given a 2x2 arrangement of frequencies (either implemented as a "table", a "matrix"
#' or a "numeric" object):
#' 
#' | n_11 | n_10|
#' |------|-----|
#' | n_01 | n_00|
#'
#' The subindex '11' corresponds to those
#' GO items enriched in both lists, '01' to items enriched in the second list but not in the first one,
#' '10' to items enriched in the first list but not enriched in the second one and '00' corresponds
#' to those GO items non enriched in both gene lists, i.e., to the double negatives, a value which
#' is ignored in the computations, except if \code{boot == TRUE}.
#'
#' In the "numeric" interface, if \code{length(x) >= 4}, the values are interpreted
#' as
#' \eqn{(n_{11}, n_{01}, n_{10}), n_{00})}{%
#' (n_11, n_01, n_10, n_00)}, always in this order and discarding extra values if necessary.
#'
#' Arguments \code{dis}, \code{se} and \code{z.conf.level} are not required. If known in advance (e.g., as
#' a consequence of previous computations with the same data), providing its value may speed the computations.
#'
#' By default, \code{z.conf.level} corresponds to the 1 - conf.level quantile of a standard normal N(0,1)
#' distribution, as the studentized statistic (^d - d) / ^se) is asymptotically N(0,1). In
#' the studentized statistic, d stands for the "true" Sorensen-Dice dissimilarity, ^d to its sample estimate
#' and ^se for the estimate of its standard error.
#' In fact, the normal is its limiting distribution but, for finite samples, the true sampling
#' distribution may present departures from normality (mainly with some inflation in the
#' left tail).
#' Computing in advance the value of argument \code{z.conf.level} may be a way to cope with
#' these departures from normality, by means of a more adequate quantile function.
#' Alternatively, if \code{boot == TRUE}, a bootstrap quantile is internally computed.
#'
#' If \code{x} is an object of class "character", it must represent a list of gene identifiers. Then the
#' confidence interval for the dissimilarity between lists \code{x} and \code{y} is computed, after summarizing
#' these gene lists as a 2x2 contingency table of joint enrichment. This last task is performed by function
#' \code{\link{buildEnrichTable}}. In general, the argument \code{...} provides extra information for this last
#' function, like p and q cut-off values to stablish enrichment.
#'
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
#' data(tab_atlas.sanger_BP3)
#' tab_atlas.sanger_BP3
#' # Gene universe:
#' data(humanEntrezIDs)
#' duppSorensen(tab_atlas.sanger_BP3)
#' dSorensen(tab_atlas.sanger_BP3) + qnorm(0.95) * seSorensen(tab_atlas.sanger_BP3)
#' # Using the bootstrap approximation instead of the normal approximation to
#' # the sampling distribution of (^d - d) / se(^d):
#' duppSorensen(tab_atlas.sanger_BP3, boot = TRUE)
#'
#' # Contingency table as a numeric vector:
#' duppSorensen(c(56, 1, 30, 47))
#' duppSorensen(c(56, 1, 30))
#'
#' # Upper confidence limit for the Sorensen-Dice dissimilarity, from scratch, 
#' # directly from two gene lists:
#' # (These examples may be considerably time consuming due to many enrichment 
#' # tests to build the contingency tables of mutual enrichment)
#' # ?pbtGeneLists
#' # duppSorensen(pbtGeneLists[[2]], pbtGeneLists[[4]],
#' #              listNames = names(pbtGeneLists)[c(2,4)],
#' #              onto = "CC", GOLevel = 5,
#' #              geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
#' # Even more time consuming (all pairwise values):
#' # duppSorensen(pbtGeneLists,
#' #              onto = "CC", GOLevel = 5,
#' #              geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
#'
#' @importFrom stats qnorm quantile rmultinom
#' @md

#' @export
duppSorensen <- function(x, ...) {
  UseMethod("duppSorensen")
}

#' @describeIn duppSorensen S3 method for class "table"
#' @export
duppSorensen.table <- function(x,
                               dis = dSorensen.table(x, check.table = FALSE),
                               se = seSorensen.table(x, check.table = FALSE),
                               conf.level = 0.95, z.conf.level = qnorm(1 - conf.level),
                               boot = FALSE, nboot = 10000,
                               listNames = NULL, check.table = TRUE, ...){
  if (check.table){
    if (!nice2x2Table(x)) {
      print(x)
      stop("Inadequate contingency table")}
  }
  if ((se == 0.0) || !is.finite(dis) || !is.finite(se)) {
    warning("Null standard error or not finite dissimilarity or standard error")
    return(NA)
  }
  if (boot) {
    n <- sum(x)
    pTab <- x / n
    bootTabs <- rmultinom(nboot, size = n, prob = pTab)
    tStats <- apply(bootTabs, 2, boot.tStat, dis = dis)
    tStats <- tStats[is.finite(tStats)]
    if (length(tStats) < nboot) {
      warMsg <- paste0("Non finite values generated in the bootstrap process (",
                       length(tStats), " out of ", nboot, ")")
      warning(warMsg)
    }
    z.conf.level <- quantile(tStats, probs = 1 - conf.level)
  }
  result <- min(dis - z.conf.level * se, 1.0)
  if (boot) {
    attr(result, "nboot") <- length(tStats)
  }
  if (is.null(listNames)) {
    names(result) <- NULL
  } else {
    names(result) <- paste("Sorensen-Dice disimilarity upper confidence limit ", listNames[1], ",", listNames[2], sep = "")
  }
  return(result)
}

#' @describeIn duppSorensen S3 method for class "matrix"
#' @export
duppSorensen.matrix <- function(x,
                                dis = dSorensen.matrix(x, check.table = FALSE),
                                se = seSorensen.matrix(x, check.table = FALSE),
                                conf.level = 0.95, z.conf.level = qnorm(1 - conf.level),
                                boot = FALSE, nboot = 10000,
                                listNames = NULL, check.table = TRUE, ...){
  if (check.table){
    if (!nice2x2Table(x)) {
      print(x)
      stop("Inadequate contingency table")}
  }
  if ((se == 0.0) || !is.finite(dis) || !is.finite(se)) {
    warning("Null standard error or not finite dissimilarity or standard error")
    return(NA)
  }
  if (boot) {
    n <- sum(x)
    pTab <- x / n
    bootTabs <- rmultinom(nboot, size = n, prob = pTab)
    tStats <- apply(bootTabs, 2, boot.tStat, dis = dis)
    tStats <- tStats[is.finite(tStats)]
    if (length(tStats) < nboot) {
      warning("Non finite values generated in the bootstrap process")
    }
    z.conf.level <- quantile(tStats, probs = 1 - conf.level)
  }
  result <- min(dis - z.conf.level * se, 1.0)
  if (boot) {
    attr(result, "nboot") <- length(tStats)
  }
  if (is.null(listNames)) {
    names(result) <- NULL
  } else {
    names(result) <- paste("Sorensen-Dice disimilarity upper confidence limit ", listNames[1], ",", listNames[2], sep = "")
  }
  return(result)
}

#' @describeIn duppSorensen S3 method for class "numeric"
#' @export
duppSorensen.numeric <- function(x,
                                 dis = dSorensen.numeric(x, check.table = FALSE),
                                 se = seSorensen.numeric(x, check.table = FALSE),
                                 conf.level = 0.95, z.conf.level = qnorm(1 - conf.level),
                                 boot = FALSE, nboot = 10000,
                                 listNames = NULL, check.table = TRUE, ...){
  if (check.table){
    if (!nice2x2Table(x)) {
      print(x)
      stop("Inadequate contingency table")}
  }
  if ((se == 0.0) || !is.finite(dis) || !is.finite(se)) {
    warning("Null standard error or not finite dissimilarity or standard error")
    return(NA)
  }
  if (boot) {
    if (length(x) < 4) {
      print(x)
      message("A numeric vector of almost 4 frequencies is required to bootstrap")
      stop()
    }
    n <- sum(x)
    pTab <- x / n
    bootTabs <- rmultinom(nboot, size = n, prob = pTab)
    tStats <- apply(bootTabs, 2, boot.tStat, dis = dis)
    tStats <- tStats[is.finite(tStats)]
    if (length(tStats) < nboot) {
      warning("Non finite values generated in the bootstrap process")
    }
    z.conf.level <- quantile(tStats, probs = 1 - conf.level)
  }
  result <- min(dis - z.conf.level * se, 1.0)
  if (boot) {
    attr(result, "nboot") <- length(tStats)
  }
  if (is.null(listNames)) {
    names(result) <- NULL
  } else {
    names(result) <- paste("Sorensen-Dice disimilarity upper confidence limit ", listNames[1], ",", listNames[2], sep = "")
  }
  return(result)
}

#' @describeIn duppSorensen S3 method for class "character"
#' @export
duppSorensen.character <- function(x, y,
                                   conf.level = 0.95,
                                   boot = FALSE, nboot = 10000,
                                   listNames = NULL, check.table = TRUE,
                                   ...){
  tab <- buildEnrichTable(x, y, listNames, check.table, ...)
  # Typical ... arguments:
  # geneUniverse=humanEntrezIDs, orgPackg="org.Hs.eg.db",
  # onto = onto, GOLevel = ontoLevel,
  result <- duppSorensen.table(tab,
                               conf.level = conf.level,
                               boot = boot, nboot = nboot,
                               listNames = listNames, check.table = check.table)
  # if (is.null(listNames)) {
  #   names(result) <- NULL
  # } else {
  #   names(result) <- paste("Sorensen-Dice disimilarity upper confidence limit ", listNames[1], ",", listNames[2], sep = "")
  # }
  return(result)
}

#' @describeIn duppSorensen S3 method for class "list"
#' @export
duppSorensen.list <- function(x,
                              conf.level = 0.95, check.table = TRUE,
                              boot = FALSE, nboot = 10000,
                              ...){
  numLists <- length(x)
  lstNams <- names(x)
  result <- matrix(0.0, ncol = numLists, nrow = numLists)
  for (iLst1 in 2:numLists) {
    for (iLst2 in 1:(iLst1-1)) {
      result[iLst1, iLst2] <- duppSorensen.character(x[[iLst1]], x[[iLst2]],
                                                     conf.level = conf.level,
                                                     boot = boot, nboot = nboot,
                                                     listNames = NULL, check.table = check.table, ...)
    }
  }
  result[upper.tri(result)] <- t(result)[upper.tri(result)]
  colnames(result) <- lstNams
  rownames(result) <- lstNams
  return(result)
}
