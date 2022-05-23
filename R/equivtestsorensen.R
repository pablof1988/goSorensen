#' Equivalence test based on the Sorensen-Dice dissimilarity
#'
#' Equivalence test based on the Sorensen-Dice dissimilarity, computed either by an asymptotic
#' normal or by a bootstrap approach.
#'
#' @param x either an object of class "table", "matrix", "numeric", "character" or "list".
#' See the details section for more information.
#' @param y an object of class "character" representing a list of gene identifiers.
#' @param d0 equivalence threshold for the Sorensen-Dice dissimilarity, d.
#' The null hypothesis states that d >= d0, i.e., inequivalence between the compared
#' gene lists and the alternative that d < d0, i.e., equivalence or dissimilarity
#' irrelevance (up to a level d0).
#' @param conf.level confidence level of the one-sided confidence interval, a value between 0 and 1.
#' @param boot boolean. If TRUE, the confidence interval and the test p-value are computed by means
#' of a bootstrap approach instead of the asymptotic normal approach. Defaults to FALSE.
#' @param nboot numeric, number of bootstrap replicates. Ignored if \code{boot == FALSE}. Defaults to 10000.
#' @param check.table Boolean. If TRUE (default), argument \code{x} is checked to adequately
#' represent a 2x2 contingency table. This checking is performed by means of function
#' \code{nice2x2Table}.
#' @param listNames character(2), names of both gene lists to be compared (only for the character interface)
#' @param ... extra parameters for function \code{buildEnrichTable}.
#'
#' @return
#' For all interfaces (except for the "list" interface) the result is a list of class "equivSDhtest" which
#' inherits from "htest", with the following components:
#' \describe{
#'   \item{statistic}{the value of the studentized statistic (dSorensen(x) - d0) / seSorensen(x)}
#'   \item{p.value}{the p-value of the test}
#'   \item{conf.int}{the one-sided confidence interval (0, dUpp]}
#'   \item{estimate}{the Sorensen dissimilarity estimate, dSorensen(x)}
#'   \item{null.value}{the value of d0}
#'   \item{stderr}{the standard error of the Sorensen dissimilarity estimate, seSorensen(x),
#'   used as denominator in the studentized statistic}
#'   \item{alternative}{a character string describing the alternative hypothesis}
#'   \item{method}{a character string describing the test}
#'   \item{data.name}{a character string giving the names of the data}
#'   \item{enrichTab}{the 2x2 contingency table of joint enrichment whereby the test was based}
#' }
#' For the "list" interface, the result is an "equivSDhtestList", a list of objects with
#' all pairwise comparisons, each one being an object of "equivSDhtest" class.
#'
#' @details
#' This function computes either the normal asymptotic or the bootstrap equivalence
#' test based on the Sorensen-Dice dissimilarity, given a 2x2 arrangement of frequencies
#' (either implemented as a "table", a "matrix" or a "numeric" object):
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
#' \eqn{(n_{11}, n_{01}, n_{10}, n_{00})}{%
#' (n_11, n_01, n_10, n_00)}, always in this order and discarding extra values if necessary.
#'
#' If \code{x} is an object of class "character", it must represent a list of gene identifiers. Then the
#' equivalence test is performed between lists \code{x} and \code{y}, after internally summarizing these
#' gene lists as a 2x2 contingency table of joint enrichment.
#'
#' If \code{x} is an object of class "list", each of its elements must be a "character" vector of gene
#' identifiers. Then all pairwise equivalence tests are performed between these gene lists.

#' The test is based on the fact that the studentized statistic (^d - d) / ^se is aproximately distributed
#' as a standard normal. ^d stands for the sample Sorensen-Dice dissimilarity, d for its true (unknown) value
#' and ^se for the estimate of its standard error. This result is asymptotically correct, but the true distribution
#' of the studentized statistic is not exactly normal for finite samples, with a heavier left tail than expected
#' under the Gaussian model. The bootstrap method provides a better approximation to this distribuion.
#'
#' @seealso \code{\link{nice2x2Table}} for checking and reformatting data,
#' \code{\link{dSorensen}} for computing the Sorensen-Dice dissimilarity,
#' \code{\link{seSorensen}} for computing the standard error of the dissimilarity,
#' \code{\link{duppSorensen}} for the upper limit of a one-sided confidence interval
#' of the dissimilarity.
#' \code{\link{getTable}}, \code{\link{getPvalue}}, \code{\link{getUpper}},
#' \code{\link{getSE}} for accessing specific fields in the result of these testing functions.
#' \code{\link{update}} for updating the result of these testing functions with alternative
#' equivalence limits or confidence levels.
#'
#' @examples
#' # Gene lists 'atlas' and 'sanger' in 'allOncoGeneLists' dataset. Table of joint enrichment
#' # of GO items in ontology BP at level 3.
#' tab_atlas.sanger_BP3
#' equivTestSorensen(tab_atlas.sanger_BP3)
#' # Bootstrap test:
#' equivTestSorensen(tab_atlas.sanger_BP3, boot = TRUE)
#'
#' # Equivalence tests from scratch, directly from gene lists:
#' # (These examples may be considerably time consuming due to many enrichment
#' # tests to build the contingency tables of mutual enrichment)
#' # ?pbtGeneLists
#' # Gene universe:
#' # data(humanEntrezIDs)
#' # equivTestSorensen(pbtGeneLists[["IRITD3"]], pbtGeneLists[["IRITD5"]],
#' #                   geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#' #                   onto = "CC", GOLevel = 5)
#' # Bootstrap instead of normal approximation test:
#' # equivTestSorensen(pbtGeneLists[["IRITD3"]], pbtGeneLists[["IRITD5"]],
#' #                   geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#' #                   onto = "CC", GOLevel = 5,
#' #                   boot = TRUE)
#'
#' # Essentially, the above code makes:
#' # IRITD3vs5.CC5 <- buildEnrichTable(pbtGeneLists[["IRITD3"]], pbtGeneLists[["IRITD5"]],
#' #                                   geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#' #                                   onto = "CC", GOLevel = 5)
#' # IRITD3vs5.CC5
#' # equivTestSorensen(IRITD3vs5.CC5)
#' # equivTestSorensen(IRITD3vs5.CC5, boot = TRUE)
#' # (Note that building first the contingency table may be advantageous to save time!)
#'
#' # All pairwise equivalence tests:
#' # equivTestSorensen(pbtGeneLists,
#' #                   geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#' #                   onto = "CC", GOLevel = 5)
#'
#'
#' # Equivalence test on a contingency table represented as a numeric vector:
#' equivTestSorensen(c(56, 1, 30, 47))
#' equivTestSorensen(c(56, 1, 30, 47), boot = TRUE)
#' equivTestSorensen(c(56, 1, 30))
#' # Error: all frequencies are needed for bootstrap:
#' try(equivTestSorensen(c(56, 1, 30), boot = TRUE), TRUE)
#' @md

#' @importFrom stats pnorm qnorm quantile rmultinom
#' @export
equivTestSorensen <- function(x, ...) {
  UseMethod("equivTestSorensen")
}

#' @describeIn equivTestSorensen S3 method for class "table"
#' @export
equivTestSorensen.table <- function(x,
                                    d0 = 1 / (1 + 1.25),
                                    conf.level = 0.95,
                                    boot = FALSE, nboot = 10000,
                                    check.table = TRUE, ...){
  data.name <- deparse(substitute(x))
  if (check.table){
    if (!nice2x2Table.table(x)) {
      print(x)
      stop("Inadequate table to compute the standard error")
    }
  }
  se <- seSorensen.table(x, check.table = FALSE)
  d <- dSorensen.table(x, check.table = FALSE)
  names(d) <- "Sorensen dissimilarity"
  names(se) <- "standard error"
  attr(d, "se") <- se
  names(d0) <- "equivalence limit d0"
  stat <- (d - d0) / se
  names(stat) <- "(d - d0) / se"
  if (!is.finite(stat)) {
    p.val <- NA
    meth <- "No test performed due not finite (d - d0) / se statistic"
    conf.int <- c(0.0, NA)
  } else {
    if (boot) {
      n <- sum(x)
      pTab <- x / n
      bootTabs <- rmultinom(nboot, size = n, prob = pTab)
      tStats <- apply(bootTabs, 2, boot.tStat, dis = d)
      tStats <- tStats[is.finite(tStats)]
      len.tStats <- length(tStats)
      if (len.tStats < nboot) {
        warning("Non finite values generated in the bootstrap process")
      }
      z.conf.level <- quantile(tStats, probs = 1 - conf.level)
      p.val <- (sum(tStats <= stat) + 1) / (len.tStats + 1)
      meth <- paste0(
        "Bootstrap test for 2x2 contingency tables based on the Sorensen-Dice dissimilarity (nboot = ",
        len.tStats, ")")
      attr(meth, "nboot") <- len.tStats
    } else {
        z.conf.level <- qnorm(1 - conf.level)
        p.val <- pnorm(stat)
        meth <- "Normal asymptotic test for 2x2 contingency tables based on the Sorensen-Dice dissimilarity"
    }
    du <- d - z.conf.level * se
    conf.int <- c(0, min(1, du))
  }
  names(p.val) <- "p-value"
  attr(conf.int, "conf.level") <- conf.level
  names(conf.int) <- c("confidence interval", "dUpper")
  result <- list(statistic = stat,
                 p.value = p.val,
                 conf.int = conf.int,
                 estimate = d,
                 null.value = d0,
                 stderr = se,
                 alternative = "less",
                 method = meth,
                 data.name = data.name,
                 enrichTab = x)
  attr(result, "check.table") <- check.table
  class(result) <- c("equivSDhtest", "htest")
  return(result)
}

#' @describeIn equivTestSorensen S3 method for class "matrix"
#' @export
equivTestSorensen.matrix <- function(x,
                                     d0 = 1 / (1 + 1.25),
                                     conf.level = 0.95,
                                     boot = FALSE, nboot = 10000,
                                     check.table = TRUE, ...){
  data.name <- deparse(substitute(x))
  if (check.table){
    if (!nice2x2Table.matrix(x)) {
      print(x)
      stop("Inadequate table to compute the standard error")
    }
  }
  se <- seSorensen.matrix(x, check.table = FALSE)
  d <- dSorensen.matrix(x, check.table = FALSE)
  names(d) <- "Sorensen dissimilarity"
  names(se) <- "standard error"
  attr(d, "se") <- se
  names(d0) <- "equivalence limit d0"
  stat <- (d - d0) / se
  names(stat) <- "(d - d0) / se"
  if (!is.finite(stat)) {
    p.val <- NA
    meth <- "No test performed due not finite (d - d0) / se statistic"
    conf.int <- c(0.0, NA)
  } else {
    if (boot) {
      n <- sum(x)
      pTab <- x / n
      bootTabs <- rmultinom(nboot, size = n, prob = pTab)
      tStats <- apply(bootTabs, 2, boot.tStat, dis = d)
      tStats <- tStats[is.finite(tStats)]
      len.tStats <- length(tStats)
      if (len.tStats < nboot) {
        warning("Non finite values generated in the bootstrap process")
      }
      z.conf.level <- quantile(tStats, probs = 1 - conf.level)
      p.val <- (sum(tStats <= stat) + 1) / (len.tStats + 1)
      meth <- paste0(
        "Bootstrap test for 2x2 contingency tables based on the Sorensen-Dice dissimilarity (nboot = ",
        len.tStats, ")")
      attr(meth, "nboot") <- len.tStats
    } else {
      z.conf.level <- qnorm(1 - conf.level)
      p.val <- pnorm(stat)
      meth <- "Normal asymptotic test for 2x2 contingency tables based on the Sorensen-Dice dissimilarity"
    }
    du <- d - z.conf.level * se
    conf.int <- c(0, min(1, du))
  }
  names(p.val) <- "p-value"
  attr(conf.int, "conf.level") <- conf.level
  names(conf.int) <- c("confidence interval", "dUpper")
  result <- list(statistic = stat,
                 p.value = p.val,
                 conf.int = conf.int,
                 estimate = d,
                 null.value = d0,
                 stderr = se,
                 alternative = "less",
                 method = meth,
                 data.name = data.name,
                 enrichTab = x)
  attr(result, "check.table") <- check.table
  class(result) <- c("equivSDhtest", "htest")
  return(result)
}

#' @describeIn equivTestSorensen S3 method for class "numeric"
#' @export
equivTestSorensen.numeric <- function(x,
                                      d0 = 1 / (1 + 1.25),
                                      conf.level = 0.95,
                                      boot = FALSE, nboot = 10000,
                                      check.table = TRUE, ...){
  data.name <- deparse(substitute(x))
  if (check.table){
    if (!nice2x2Table.numeric(x)) {
      print(x)
      stop("Inadequate table to compute the standard error")
    }
  }
  se <- seSorensen.matrix(x, check.table = FALSE)
  d <- dSorensen.matrix(x, check.table = FALSE)
  names(d) <- "Sorensen dissimilarity"
  names(se) <- "standard error"
  attr(d, "se") <- se
  names(d0) <- "equivalence limit d0"
  stat <- (d - d0) / se
  names(stat) <- "(d - d0) / se"
  if (!is.finite(stat)) {
    p.val <- NA
    meth <- "No test performed due non finite (d - d0) / se statistic"
    conf.int <- c(0.0, NA)
  } else {
    if (boot) {
      if (length(x) < 4) {
        print(x)
        message("A numeric vector of almost 4 frequencies is required to bootstrap")
        stop()
      }
      n <- sum(x)
      pTab <- x / n
      bootTabs <- rmultinom(nboot, size = n, prob = pTab)
      tStats <- apply(bootTabs, 2, boot.tStat, dis = d)
      tStats <- tStats[is.finite(tStats)]
      len.tStats <- length(tStats)
      if (len.tStats < nboot) {
        warning("Non finite values generated in the bootstrap process")
      }
      z.conf.level <- quantile(tStats, probs = 1 - conf.level)
      p.val <- (sum(tStats <= stat) + 1) / (len.tStats + 1)
      meth <- paste0(
        "Bootstrap test for 2x2 contingency tables based on the Sorensen-Dice dissimilarity (nboot = ",
        len.tStats, ")")
      attr(meth, "nboot") <- len.tStats
    } else {
      z.conf.level <- qnorm(1 - conf.level)
      p.val <- pnorm(stat)
      meth <- "Normal asymptotic test for 2x2 contingency tables based on the Sorensen-Dice dissimilarity"
    }
    du <- d - z.conf.level * se
    conf.int <- c(0, min(1, du))
  }
  names(p.val) <- "p-value"
  attr(conf.int, "conf.level") <- conf.level
  names(conf.int) <- c("confidence interval", "dUpper")
  result <- list(statistic = stat,
                 p.value = p.val,
                 conf.int = conf.int,
                 estimate = d,
                 null.value = d0,
                 stderr = se,
                 alternative = "less",
                 method = meth,
                 data.name = data.name,
                 enrichTab = x)
  attr(result, "check.table") <- check.table
  class(result) <- c("equivSDhtest", "htest")
  return(result)

}

#' @describeIn equivTestSorensen S3 method for class "character"
#' @export
equivTestSorensen.character <- function(x, y, d0 = 1 / (1 + 1.25),
                                        conf.level = 0.95,
                                        boot = FALSE, nboot = 10000,
                                        listNames = c("gene.list1", "gene.list2"),
                                        check.table = TRUE,
                                        ...){
  tab <- buildEnrichTable(x, y, listNames, check.table = FALSE, ...)
  # Typical ... arguments:
  # geneUniverse=humanEntrezIDs, orgPackg="org.Hs.eg.db",
  # onto = onto, GOLevel = ontoLevel,
  return(equivTestSorensen.table(tab, d0 = d0,
                                 boot = boot, nboot = nboot,
                                 conf.level = conf.level, check.table = check.table))
}

#' @describeIn equivTestSorensen S3 method for class "list"
#' @export
equivTestSorensen.list <- function(x, d0 = 1 / (1 + 1.25),
                                   conf.level = 0.95, boot = FALSE, nboot = 10000,
                                   check.table = TRUE,
                                   ...){
  numLists <- length(x)
  lstNams <- names(x)
  equivTests <- lapply(2:numLists, function(iLst1) {
    oneVsOthers <- lapply(1:(iLst1-1), function(iLst2) {
      return(equivTestSorensen.character(x[[iLst1]], x[[iLst2]],
                                         d0 = d0, conf.level = conf.level,
                                         boot = boot, nboot = nboot,
                                         check.table = check.table,
                                         listNames = c(lstNams[iLst1], lstNams[iLst2]),
                                         ...))
    })
    names(oneVsOthers) <- lstNams[1:(iLst1-1)]
    return(oneVsOthers)
  })
  names(equivTests) <- lstNams[2:numLists]
  class(equivTests) <- c("equivSDhtestList", "list")
  return(equivTests)
}

#' Iterate \code{equivTestSorensen} along the specified GO ontologies and GO levels
#'
#' @param x object of class "list". Each of its elements must be a "character" vector of gene
#' identifiers. Then all pairwise equivalence tests are performed between these gene lists,
#' iterating the process for all specified GO ontologies and GO levels.
#' @param d0 equivalence threshold for the population Sorensen-Dice dissimilarity, d. The null hypothesis
#' states that d >= d0, and the alternative that d < d0.
#' @param conf.level confidence level of the one-sided confidence interval, a value between 0 and 1.
#' @param boot boolean. If TRUE, the confidence interval and the test p-value are computed by means
#' of a bootstrap approach instead of the asymptotic normal approach. Defaults to FALSE.
#' @param nboot numeric, number of bootstrap replicates. Ignored if \code{boot == FALSE}. Defaults to 10000.
#' @param check.table Boolean. If TRUE (default), argument \code{x} is checked to adequately
#' represent a 2x2 contingency table. This checking is performed by means of function
#' \code{nice2x2Table}.
#' @param ontos "character", GO ontologies to analyse. Defaults to \code{c("BP", "CC", "MF")}.
#' @param GOLevels "integer", GO levels to analyse inside each of these GO ontologies.
#' @param ... extra parameters for function \code{buildEnrichTable}.
#'
#' @return
#' An object of class "AllEquivSDhtest". It is a list with as many components as GO ontologies have been analysed.
#' Each of these elements is itself a list with as many components as GO levels have been analized.
#' Finally, the elements of these lists are objects as generated by \code{equivTestSorensen.list},
#' i.e., objects of class "equivSDhtestList" containing all pairwise comparisons between the gene lists
#' in argument \code{x}.
#'
#' @examples
#' # This example is extremely time consuming, it scans two GO ontologies and three
#' # GO levels inside them to perform the equivalence test.
#' # Gene universe:
#' # data(humanEntrezIDs)
#' # Gene lists to be explored for enrichment:
#' # data(pbtGeneLists)
#' # allEquivTestSorensen(pbtGeneLists,
#' #                      geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#' #                      ontos = c("MF", "BP"), GOLevels = 4:6)
#'
#' @export
allEquivTestSorensen <- function(x, d0 = 1 / (1 + 1.25), conf.level = 0.95,
                                 boot = FALSE, nboot = 10000, check.table = TRUE,
                                 ontos = c("BP", "CC", "MF"),
                                 GOLevels = 3:10,
                                 ...){
  allOntos <- lapply(ontos, function(onto) {
    thisOnto <- lapply(GOLevels, function(lev) {
      cat("========================================================\n")
      cat("Ontology ", onto, " level ", lev, "\n")
      result <- equivTestSorensen(x, d0 = d0, conf.level = conf.level,
                                  boot = boot, nboot = nboot,
                                  check.table = check.table,
                                  onto = onto, GOLevel = lev, ...)
      print(result)
      return(result)
    })
    names(thisOnto) <- paste("level", GOLevels)
    return(thisOnto)
  })
  names(allOntos) <- ontos
  class(allOntos) <- c("AllEquivSDhtest", "list")
  return(allOntos)
}
