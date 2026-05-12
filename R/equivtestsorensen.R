#' Equivalence test based on the Sorensen-Dice dissimilarity
#'
#' Equivalence test based on the Sorensen-Dice dissimilarity, computed either by
#' an asymptotic normal approach or by a bootstrap approach.
#'
#' @param x either an object of class "table", "matrix", "numeric", "character",
#' "list" or "tableList".
#' See the details section for more information.
#' @return See method-specific documentation.
#' @export
equivTestSorensen <- function(x, ...) {
  UseMethod("equivTestSorensen")
}

#' @export
equivTestSorensen.table <- function(x,
                                    d0 = 1 / (1 + 1.25),
                                    conf.level = 0.95,
                                    boot = FALSE, nboot = 10000,
                                    check.table = TRUE, ...) {
  stopifnot(
    "Arguments 'd0' and 'conf.level' must be numeric between 0 and 1" =
      (0 < d0) && (d0 < 1) && (0 < conf.level) && (conf.level < 1)
  )
  stopifnot("Argument 'check.table' must be logical" = is.logical(check.table))
  if (check.table) {
    nice2x2Table.table(x)
  }
  data.name <- deparse(substitute(x))
  se <- seSorensen.table(x, check.table = FALSE)
  d <- dSorensen.table(x, check.table = FALSE)
  names(d) <- "Sorensen dissimilarity"
  names(se) <- "standard error"
  attr(d, "se") <- se
  names(d0) <- "equivalence limit d0"
  stat <- (d - d0) / se
  names(stat) <- "(d - d0) / se"
  if (!is.finite(stat)) {
    p.val <- NaN
    meth <- "No test performed due not finite (d - d0) / se statistic"
    conf.int <- c(0.0, NaN)
    if (boot) {
      len.tStats <- NaN
    }
  } else {
    if (boot) {
      stopifnot("Argument 'nboot' must be numeric" = is.numeric(nboot))
      n <- sum(x)
      pTab <- x / n
      bootTabs <- rmultinom(nboot, size = n, prob = pTab)
      tStats <- apply(bootTabs, 2, boot.tStat, dis = d)
      tStats <- tStats[is.finite(tStats)]
      len.tStats <- length(tStats)
      z.conf.level <- quantile(tStats, probs = 1 - conf.level)
      p.val <- (sum(tStats <= stat) + 1) / (len.tStats + 1)
      meth <- "Bootstrap test for 2x2 contingency tables based on the Sorensen-
      Dice dissimilarity"
      if (len.tStats < nboot) {
        meth <- paste(meth, "\n(", len.tStats, " effective bootstrap replicates
                      of ", nboot, ")", sep = "")
      } else {
        meth <- paste(meth, "\n(", nboot, " bootstrap replicates)", sep = "")
      }
    } else {
      z.conf.level <- qnorm(1 - conf.level)
      p.val <- pnorm(stat)
      meth <- "Normal asymptotic test for 2x2 contingency tables based on the
        Sorensen-Dice dissimilarity"
    }
    du <- d - z.conf.level * se
    conf.int <- c(0, min(1, du))
  }
  names(p.val) <- "p-value"
  attr(conf.int, "conf.level") <- conf.level
  names(conf.int) <- c("confidence interval", "dUpper")
  result <- list(
    statistic = stat,
    p.value = p.val,
    conf.int = conf.int,
    estimate = d,
    null.value = d0,
    stderr = se,
    alternative = "less",
    method = meth,
    data.name = data.name,
    enrichTab = x
  )
  attr(result, "check.table") <- check.table
  if (boot) {
    attr(result, "nboot") <- nboot
    attr(result, "eff.nboot") <- len.tStats
  }
  class(result) <- c("equivSDhtest", "htest")
  return(result)
}


#' @export
equivTestSorensen.matrix <- function(x,
                                     d0 = 1 / (1 + 1.25),
                                     conf.level = 0.95,
                                     boot = FALSE, nboot = 10000,
                                     check.table = TRUE, ...) {
  stopifnot(
    "Arguments 'd0' and 'conf.level' must be numeric between 0 and 1" =
      (0 < d0) && (d0 < 1) && (0 < conf.level) && (conf.level < 1)
  )
  stopifnot("Argument 'check.table' must be logical" = is.logical(check.table))
  if (check.table) {
    nice2x2Table.matrix(x)
  }
  data.name <- deparse(substitute(x))
  se <- seSorensen.matrix(x, check.table = FALSE)
  d <- dSorensen.matrix(x, check.table = FALSE)
  names(d) <- "Sorensen dissimilarity"
  names(se) <- "standard error"
  attr(d, "se") <- se
  names(d0) <- "equivalence limit d0"
  stat <- (d - d0) / se
  names(stat) <- "(d - d0) / se"
  if (!is.finite(stat)) {
    p.val <- NaN
    meth <- "No test performed due not finite (d - d0) / se statistic"
    conf.int <- c(0.0, NaN)
    if (boot) {
      len.tStats <- NaN
    }
  } else {
    if (boot) {
      stopifnot("Argument 'nboot' must be numeric" = is.numeric(nboot))
      n <- sum(x)
      pTab <- x / n
      bootTabs <- rmultinom(nboot, size = n, prob = pTab)
      tStats <- apply(bootTabs, 2, boot.tStat, dis = d)
      tStats <- tStats[is.finite(tStats)]
      len.tStats <- length(tStats)
      z.conf.level <- quantile(tStats, probs = 1 - conf.level)
      p.val <- (sum(tStats <= stat) + 1) / (len.tStats + 1)
      meth <- "Bootstrap test for 2x2 contingency tables based on the Sorensen-
      Dice dissimilarity"
      if (len.tStats < nboot) {
        meth <- paste(meth, "\n(", len.tStats, " effective bootstrap replicates
                      of ", nboot, ")", sep = "")
      } else {
        meth <- paste(meth, "\n(", nboot, " bootstrap replicates)", sep = "")
      }
    } else {
      z.conf.level <- qnorm(1 - conf.level)
      p.val <- pnorm(stat)
      meth <- "Normal asymptotic test for 2x2 contingency tables based on the
      Sorensen-Dice dissimilarity"
    }
    du <- d - z.conf.level * se
    conf.int <- c(0, min(1, du))
  }
  names(p.val) <- "p-value"
  attr(conf.int, "conf.level") <- conf.level
  names(conf.int) <- c("confidence interval", "dUpper")
  result <- list(
    statistic = stat,
    p.value = p.val,
    conf.int = conf.int,
    estimate = d,
    null.value = d0,
    stderr = se,
    alternative = "less",
    method = meth,
    data.name = data.name,
    enrichTab = x
  )
  attr(result, "check.table") <- check.table
  if (boot) {
    attr(result, "nboot") <- nboot
    attr(result, "eff.nboot") <- len.tStats
  }
  class(result) <- c("equivSDhtest", "htest")
  return(result)
}


#' @export
equivTestSorensen.numeric <- function(x,
                                      d0 = 1 / (1 + 1.25),
                                      conf.level = 0.95,
                                      boot = FALSE, nboot = 10000,
                                      check.table = TRUE, ...) {
  stopifnot(
    "Arguments 'd0' and 'conf.level' must be numeric between 0 and 1" =
      (0 < d0) && (d0 < 1) && (0 < conf.level) && (conf.level < 1)
  )
  stopifnot("Argument 'check.table' must be logical" = is.logical(check.table))
  if (check.table) {
    nice2x2Table.numeric(x)
  }
  data.name <- deparse(substitute(x))
  se <- seSorensen.matrix(x, check.table = FALSE)
  d <- dSorensen.matrix(x, check.table = FALSE)
  names(d) <- "Sorensen dissimilarity"
  names(se) <- "standard error"
  attr(d, "se") <- se
  names(d0) <- "equivalence limit d0"
  stat <- (d - d0) / se
  names(stat) <- "(d - d0) / se"
  if (!is.finite(stat)) {
    p.val <- NaN
    meth <- "No test performed due non finite (d - d0) / se statistic"
    conf.int <- c(0.0, NaN)
    if (boot) {
      len.tStats <- NaN
    }
  } else {
    if (boot) {
      stopifnot(
        "Bootstraping requires a numeric vector of 4 frequencies" =
          length(x) == 4
      )
      stopifnot("Argument 'nboot' must be numeric" = is.numeric(nboot))
      n <- sum(x)
      pTab <- x / n
      bootTabs <- rmultinom(nboot, size = n, prob = pTab)
      tStats <- apply(bootTabs, 2, boot.tStat, dis = d)
      tStats <- tStats[is.finite(tStats)]
      len.tStats <- length(tStats)
      z.conf.level <- quantile(tStats, probs = 1 - conf.level)
      p.val <- (sum(tStats <= stat) + 1) / (len.tStats + 1)
      meth <- "Bootstrap test for 2x2 contingency tables based on the Sorensen-
      Dice dissimilarity"
      if (len.tStats < nboot) {
        meth <- paste(meth, "\n(", len.tStats, " effective bootstrap replicates
                      of ", nboot, ")", sep = "")
      } else {
        meth <- paste(meth, "\n(", nboot, " bootstrap replicates)", sep = "")
      }
    } else {
      z.conf.level <- qnorm(1 - conf.level)
      p.val <- pnorm(stat)
      meth <- "Normal asymptotic test for 2x2 contingency tables based on the
      Sorensen-Dice dissimilarity"
    }
    du <- d - z.conf.level * se
    conf.int <- c(0, min(1, du))
  }
  names(p.val) <- "p-value"
  attr(conf.int, "conf.level") <- conf.level
  names(conf.int) <- c("confidence interval", "dUpper")
  result <- list(
    statistic = stat,
    p.value = p.val,
    conf.int = conf.int,
    estimate = d,
    null.value = d0,
    stderr = se,
    alternative = "less",
    method = meth,
    data.name = data.name,
    enrichTab = x
  )
  attr(result, "check.table") <- check.table
  if (boot) {
    attr(result, "nboot") <- nboot
    attr(result, "eff.nboot") <- len.tStats
  }
  class(result) <- c("equivSDhtest", "htest")
  return(result)
}
#' @describeIn equivTestSorensen S3 default method.
#' @param y an object of class "character" representing a list of gene
#' identifiers (e.g., ENTREZ).
#' @param d0 equivalence threshold for the Sorensen-Dice dissimilarity, d.
#' The null hypothesis states that d >= d0, i.e., inequivalence between the
#' compared gene lists and the alternative that d < d0, i.e., equivalence or
#' dissimilarity irrelevance (up to a level d0).
#' @param conf.level confidence level of the one-sided confidence interval, a
#' value between 0 and 1.
#' @param boot boolean. If TRUE, the confidence interval and the test p-value
#' are computed by means of a bootstrap approach instead of the asymptotic
#' normal approach. Defaults to FALSE.
#' @param nboot numeric, number of initially planned bootstrap replicates.
#' Ignored if \code{boot == FALSE}. Defaults to 10000.
#' @param check.table Boolean. If TRUE (default), argument \code{x} is checked
#' to adequately represent a 2x2 contingency table (or an aggregate of them) or
#' gene lists producing a correct table. This checking is performed by means of
#' function \code{nice2x2Table}.
#' @param ... extra parameters for function \code{buildEnrichTable}.
#'
#' @return
#' For all interfaces (except for the "list" and "tableList" interfaces), the
#' result is a list of class "equivSDhtest" which inherits from "htest", with
#' the following components:
#' \describe{
#'   \item{statistic}{The value of the studentized statistic
#'   \eqn{(dSorensen(x) - d0) / seSorensen(x)}.}
#'   \item{p.value}{The p-value of the test.}
#'   \item{conf.int}{The one-sided confidence interval \eqn{(0, dUpp]}.}
#'   \item{estimate}{The Sorensen dissimilarity estimate, \eqn{dSorensen(x)}.}
#'   \item{null.value}{The value of \code{d0}.}
#'   \item{stderr}{The standard error of the Sorensen dissimilarity estimate,
#'   \eqn{seSorensen(x)}, used as denominator in the studentized statistic.}
#'   \item{alternative}{A character string describing the alternative
#'   hypothesis.}
#'   \item{method}{A character string describing the test.}
#'   \item{data.name}{A character string giving the names of the data.}
#'   \item{enrichTab}{The 2x2 contingency table of joint enrichment on which the
#'   test was based.}
#' }
#' For the "list" and "tableList" interfaces, the result is an object of class
#' "equivSDhtestList", a list of all pairwise comparisons, each one being an
#' object of class
#' "equivSDhtest".
#'
#' @details
#' This function computes either the normal asymptotic or the bootstrap
#' equivalence test based on the Sorensen-Dice dissimilarity, given a 2x2
#' arrangement of frequencies (either implemented as a "table", a "matrix" or a
#' "numeric" object):
#'
#' \tabular{rr}{
#' \eqn{n_{11}} \tab \eqn{n_{10}} \cr
#' \eqn{n_{01}} \tab \eqn{n_{00}},
#' }
#'
#' The subindex '11' corresponds to those GO terms enriched in both lists, '01'
#' to terms enriched in the second list but not in the first one, '10' to terms
#' enriched in the first list but not enriched in the second one and '00'
#' corresponds to those GO terms non enriched in both gene lists, i.e., to the
#' double negatives, a value which is ignored in the computations.
#'
#' In the "numeric" interface, if \code{length(x) >= 4}, the values are
#' interpreted as
#' \eqn{(n_{11}, n_{01}, n_{10}, n_{00})}{%
#' (n_11, n_01, n_10, n_00)}, always in this order and discarding extra values
#' if necessary.
#'
#' If \code{x} is an object of class "character", then \code{x} (and \code{y})
#' must represent two "character" vectors of valid gene identifiers
#' (e.g., ENTREZ).
#' Then the equivalence test is performed between \code{x} and \code{y}, after
#' internally summarizing them as a 2x2 contingency table of joint enrichment.
#' This last operation is performed by function \code{\link{buildEnrichTable}}
#' and "valid gene identifiers (e.g., ENTREZ)" stands for the coherency of these
#' gene identifiers with the arguments \code{geneUniverse} and \code{orgPackg}
#' of \code{buildEnrichTable}, passed by the ellipsis argument \code{...} in
#' \code{equivTestSorensen}.
#'
#' If \code{x} is an object of class "list", each of its elements must be a
#' "character" vector of gene identifiers (e.g., ENTREZ). Then all pairwise
#' equivalence tests are performed between these gene lists.
#'
#' Class "tableList" corresponds to objects representing all mutual enrichment
#' contingency tables generated in a pairwise fashion:
#' Given gene lists l1, l2, ..., lk, an object of class "tableList" (typically
#' constructed by a call to function \code{\link{buildEnrichTable}}) is a list
#' of lists of contingency tables tij generated from each pair of gene lists i
#' and j, with the following structure:
#'
#' $l2
#'
#' $l2$l1$t21
#'
#' $l3
#'
#' $l3$l1$t31, $l3$l2$t32
#'
#' ...
#'
#' $lk$l1$tk1, $lk$l2$tk2, ..., $lk$l(k-1)tk(k-1)
#'
#' If \code{x} is an object of class "tableList", the test is performed over
#' each one of these tables.
#'

#' The test is based on the fact that the studentized statistic (^d - d) / ^se
#' is approximately distributed as a standard normal. ^d stands for the sample
#' Sorensen-Dice dissimilarity, d for its true (unknown) value and ^se for the
#' estimate of its standard error.
#' This result is asymptotically correct, but the true distribution of the
#' studentized statistic is not exactly normal for finite samples, with a
#' heavier left tail than expected under the Gaussian model, which may produce
#' some type I error inflation.
#' The bootstrap method provides a better approximation to this distribution.
#' In the bootstrap approach, \code{nboot} new bootstrap contingency tables are
#' generated from a multinomial distribution with parameters
#' \code{size =} \eqn{ n = (n_{11} + n_{01} + n_{10} + n_{00})}{%
#' (n11 + n01 + n10 + n00)} and probabilities
#' \eqn{(n_{11} / n, n_{01} / n, n_{10} / n, n_{00} / n)}{%}. Sometimes, some of
#' these generated tables may present so low frequencies of enrichment that make
#' them unable for Sorensen-Dice computations.
#' As a consequence, the number of effective bootstrap samples may be lower than
#' the number of initially planned ones, \code{nboot}, but our simulation
#' studies concluded that this makes the test more conservative, less prone to
#' reject a truly false null hypothesis of inequivalence, but in any case
#' protects from inflating the type I error.
#'
#' In a bootstrap test result, use \code{getNboot} to access the
#' number of initially planned bootstrap replicates and \code{getEffNboot} to
#' access the number of finally effective bootstrap replicates.
#'
#' @seealso \code{\link{nice2x2Table}} for checking and reformatting data,
#' \code{\link{dSorensen}} for computing the Sorensen-Dice dissimilarity,
#' \code{\link{seSorensen}} for computing the standard error of the
#' dissimilarity, \code{\link{duppSorensen}} for the upper limit of a one-sided
#' confidence interval of the dissimilarity.
#' \code{\link{getTable}}, \code{\link{getPvalue}}, \code{\link{getUpper}},
#' \code{\link{getSE}}, \code{\link{getNboot}} and \code{\link{getEffNboot}} for
#' accessing specific fields in the result of these testing functions.
#' \code{\link{update}} for updating the result of these testing functions with
#' alternative equivalence limits, confidence levels or to convert a normal
#' result in a bootstrap result or the reverse.
#'
#' @examples
#' ## The following example is highly time-consuming and is therefore not run
#' ## automatically during R CMD check.
#'
#' \dontrun{
#' ## i) Obtaining ENTREZ identifiers for the gene universe of humans:
#' library(org.Hs.eg.db)
#' humanEntrezIDs <- keys(org.Hs.eg.db, keytype = "ENTREZID")
#'
#' ## ii) Gene lists to be explored for analysis:
#' data(allOncoGeneLists)
#'
#' # iii) Calculation of Calculation of the equivalence test of all joint
#' # enrichment contingency tables obtained from the BP ontology at the GO 4
#' # level.
#' eqTest_all_BP4 <- equivTestSorensen(allOncoGeneLists,
#'   geneUniverse = humanEntrezIDs,
#'   orgPackg = "org.Hs.eg.db",
#'   onto = "BP",
#'   GOLevel = 4,
#'   d0 = 0.4444,
#'   conf.level = 0.95
#' )
#' eqTest_all_BP4
#' }
#'
#' # Since running this example may take several minutes, the result has been
#' # pre-computed and is accessible as the following:
#' data(eqTest_all_BP4)
#' eqTest_all_BP4
#' # This shortcut applies only to this example; for your own gene-list data,
#' # the computation must be performed explicitly.
#'
#' # For a complete overview of this function's use, see the section 4 of the
#' # vignette "Introduction to goSorensen". You can do this by consulting the
#' # general package documentation or by directly running the following code in
#' # the R console:
#' # vignette("goSorensen_Introduction", package = "goSorensen")
#'
#' @md
#'
#' @importFrom stats pnorm qnorm quantile rmultinom
#' @export
equivTestSorensen.character <- function(x, y, d0 = 1 / (1 + 1.25),
                                        conf.level = 0.95,
                                        boot = FALSE, nboot = 10000,
                                        check.table = TRUE,
                                        ...) {
  tab <- buildEnrichTable(x, y, check.table = check.table, ...)
  # Typical ... arguments:
  # geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
  # onto = onto, GOLevel = ontoLevel,
  return(equivTestSorensen.table(tab,
    d0 = d0,
    boot = boot, nboot = nboot,
    conf.level = conf.level, check.table = FALSE
  ))
}


#' @export
equivTestSorensen.list <- function(x, d0 = 1 / (1 + 1.25),
                                   conf.level = 0.95, boot = FALSE,
                                   nboot = 10000,
                                   check.table = TRUE,
                                   ...) {
  # numLists <- length(x)
  # lstNams <- names(x)
  # equivTests <- lapply(seq.int(2, numLists), function(iLst1) {
  #   oneVsOthers <- lapply(seq.int(1, iLst1-1), function(iLst2) {
  #     return(equivTestSorensen.character(x[[iLst1]], x[[iLst2]],
  #                                        d0 = d0, conf.level = conf.level,
  #                                        boot = boot, nboot = nboot,
  #                                        check.table = check.table,
  #                                        listNames = c(lstNams[iLst1],
  #                                        lstNams[iLst2]),
  #                                        ...))
  #   })
  #   names(oneVsOthers) <- lstNams[seq.int(1, iLst1-1)]
  #   return(oneVsOthers)
  # })
  # names(equivTests) <- lstNams[seq.int(2, numLists)]
  # class(equivTests) <- c("equivSDhtestList", "list")
  # return(equivTests)
  tabs <- buildEnrichTable(x, check.table = check.table, ...)
  return(equivTestSorensen.tableList(tabs,
    d0 = d0,
    conf.level = conf.level,
    boot = boot, nboot = nboot,
    check.table = FALSE
  ))
}


#' @export
equivTestSorensen.tableList <- function(x, d0 = 1 / (1 + 1.25),
                                        conf.level = 0.95, boot = FALSE,
                                        nboot = 10000,
                                        check.table = TRUE, ...) {
  lstNams <- names(x)
  equivTests <- lapply(x, function(lst1) {
    oneVsOthers <- lapply(lst1, function(tab) {
      return(equivTestSorensen.table(tab,
        d0 = d0, conf.level = conf.level,
        boot = boot, nboot = nboot,
        check.table = check.table
      ))
    })
    names(oneVsOthers) <- names(lst1)
    return(oneVsOthers)
  })
  names(equivTests) <- lstNams
  class(equivTests) <- c("equivSDhtestList", "list")
  return(equivTests)
}
