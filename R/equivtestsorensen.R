#' Asymptotic equivalence test for the Sorensen-Dice dissimilarity
#'
#' @param x either an object of class "table", "matrix", "numeric", "character" or "list".
#' See the details section for more information.
#' @param y an object of class "character" representing a list of gene identifiers.
#' @param n total number of cases (GO nodes) in the contingency table. Only required
#' (sometimes) on the "numeric" interface, see the details section of \code{\link{nice2x2Table}}.
#' @param d0 equivalence threshold for the population Sorensen-Dice dissimilarity, d. The null hypothesis
#' states that d >= d0, and the alternative that d < d0.
#' @param conf.level confidence level of the one-sided confidence interval, a value between 0 and 1.
#' @param check.table Boolean. If TRUE (default), argument \code{x} is checked to adequately
#' represent a 2x2 contingency table. This checking is performed by means of function
#' \code{nice2x2Table}.
#' @param ... extra parameters for function \code{crossTabGOIDs4GeneLists} in package \code{equivStandardTest}.
#'
#' @return
#' For all interfaces (except for the "list" interface) the result is a list of class "equivSDhtest" descending
#' from "htest", with the following components:
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
#' In the "table", "matrix" and "numeric" interfaces, \code{x} must represent a 2x2 contingency
#' table of joint enrichment frequencies among n GO items: n00 items non enriched in both lists, n01 items
#' non enriched in the first list but enriched in the second one and so on.
#' If \code{x} is an object of class "character", it must represent a list of gene identifiers. Then the
#' equivalence test is performed between lists \code{x} and \code{y}, after internally summarizing these
#' gene lists as a 2x2 contingency table of joint enrichment.
#' If \code{x} is an object of class "list", each of its elements must be a "character" vector of gene
#' identifiers. Then all pairwise equivalence tests are performed between these gene lists.

#' The test is based on the fact that the studentized statistic is asymptotically distributed
#' as a standard normal.
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

#'
#' library(equivStandardTest)
#' # Building enrichment contingency tables from scratch, using package "equivStandardTest"
#' # Gene universe:
#' data(humanEntrezIDs)
#' # Gene lists to be explored for enrichment:
#' pbtGeneLists
#' IRITD3VS5.MF4 <- crossTabGOIDs4GeneLists(pbtGeneLists[["IRITD3"]], pbtGeneLists[["IRITD5"]],
#'                                          geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#'                                          onto = "MF", GOLevel = 4)
#' IRITD3VS5.MF4
#' equivTestSorensen(IRITD3VS5.MF4)
#' # Perform these computations in a single step:
#' equivTestSorensen(pbtGeneLists[["IRITD3"]], pbtGeneLists[["IRITD5"]],
#'                   geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#'                   onto = "MF", GOLevel = 4)
#' # All pairwise equivalence tests:
#' equivTestSorensen(pbtGeneLists,
#'                   geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#'                   onto = "MF", GOLevel = 4)

#'
#'
#' @export
equivTestSorensen <- function(x, ...) {
  UseMethod("equivTestSorensen")
}

#' @describeIn equivTestSorensen S3 method for class "table"
#' @export
equivTestSorensen.table <- function(x, n = sum(x), d0 = 1 / (1 + 1.25),
                                    conf.level = 0.95,
                                    check.table = TRUE){
  data.name <- deparse(substitute(x))
  if (check.table){
    x <- nice2x2Table(x)
  }
  d <- dSorensen.table(x)
  names(d) <- "Sorensen dissimilarity"
  names(d0) <- "equivalence limit d0"
  se <- seSorensen.table(x, n, FALSE)
  attr(d, "se") <- se
  names(se) <- "standard deviation"
  stat <- (d - d0) / se
  names(stat) <- "(d - d0) / se"
  p.val <- pnorm(stat)
  names(p.val) <- "p-value"
  du <- duppSorensen.table(x, n, d, se, conf.level, check.table = FALSE)
  conf.int <- c(0, min(1, du))
  attr(conf.int, "conf.level") <- conf.level
  names(conf.int) <- c("confidence interval", "dUpper")
  result <- list(statistic = stat,
                 p.value = p.val,
                 conf.int = conf.int,
                 estimate = d,
                 null.value = d0,
                 stderr = se,
                 alternative = "less",
                 method = "Asymptotic test for 2x2 contingency tables based on the Sorensen-Dice dissimilarity",
                 data.name = data.name,
                 enrichTab = x)
  class(result) <- c("equivSDhtest", "htest")
  return(result)
}

#' @describeIn equivTestSorensen S3 method for class "matrix"
#' @export
equivTestSorensen.matrix <- function(x, n = sum(x), d0 = 1 / (1 + 1.25),
                                     conf.level = 0.95){
  data.name <- deparse(substitute(x))
  return(equivTestSorensen.table(nice2x2Table.matrix(x),
                                 d0 = d0, conf.level = conf.level,
                                 check.table = FALSE))
  # d <- dSorensen.matrix(x)
  # names(d) <- "Sorensen dissimilarity"
  # names(d0) <- "equivalence limit d0"
  # se <- seSorensen.matrix(x, n, FALSE)
  # attr(d, "se") <- se
  # names(se) <- "standard deviation"
  # stat <- (d - d0) / se
  # names(stat) <- "(d - d0) / se"
  # p.val <- pnorm(stat)
  # names(p.val) <- "p-value"
  # du <- duppSorensen.matrix(x, n, d, se, conf.level, check.table = FALSE)
  # conf.int <- c(0, min(1, du))
  # attr(conf.int, "conf.level") <- conf.level
  # names(conf.int) <- c("confidence interval", "dUpper")
  # result <- list(statistic = stat, p.value = p.val,
  #                conf.int = conf.int,
  #                # estimate = c(d, se),
  #                estimate = d,
  #                null.value = d0, stderr = se,
  #                alternative = "less",
  #                method = "Asymptotic test for 2x2 contingency tables based on the Sorensen-Dice dissimilarity",
  #                data.name = data.name,
  #                enrichTab = x)
  # class(result) <- c("equivSDhtest", "htest")
  # return(result)
}

#' @describeIn equivTestSorensen S3 method for class "numeric"
#' @export
equivTestSorensen.numeric <- function(x, n = sum(x), d0 = 1 / (1 + 1.25),
                                      conf.level = 0.95){
  data.name <- deparse(substitute(x))
  return(equivTestSorensen.table(nice2x2Table.numeric(x, n),
                                 d0 = d0, conf.level = conf.level,
                                 check.table = FALSE))
  # d <- dSorensen.numeric(x, FALSE)
  # names(d) <- "Sorensen dissimilarity"
  # names(d0) <- "equivalence limit d0"
  # se <- seSorensen.numeric(x, n, FALSE)
  # attr(d, "se") <- se
  # names(se) <- "standard deviation"
  # stat <- (d - d0) / se
  # names(stat) <- "(d - d0) / se"
  # p.val <- pnorm(stat)
  # names(p.val) <- "p-value"
  # du <- duppSorensen.numeric(x, n, d, se, conf.level = conf.level, check.table = FALSE)
  # conf.int <- c(0, min(1, du))
  # attr(conf.int, "conf.level") <- conf.level
  # names(conf.int) <- c("confidence interval", "dUpper")
  # result <- list(statistic = stat, p.value = p.val,
  #                conf.int = conf.int,
  #                # estimate = c(d, se),
  #                estimate = d,
  #                null.value = d0, stderr = se,
  #                alternative = "less",
  #                method = "Asymptotic test for 2x2 contingency tables based on the Sorensen-Dice dissimilarity",
  #                data.name = data.name,
  #                enrichTab = x)
  # class(result) <- c("equivSDhtest", "htest")
  # return(result)
}

#' @describeIn equivTestSorensen S3 method for class "character"
#' @export
equivTestSorensen.character <- function(x, y, d0 = 1 / (1 + 1.25),
                                        conf.level = 0.95,
                                        listNames = c("gene.list1", "gene.list2"),
                                        ...){
  tab <- crossTabGOIDs4GeneLists (genelist1 = x, genelist2 = y,
                                  ...)
                                  # Typical ... arguments:
                                  # geneUniverse=humanEntrezIDs, orgPackg="org.Hs.eg.db",
                                  # onto = onto, GOLevel = ontoLevel,
  tab <- nice2x2Table.table(tab, listNames)
  return(equivTestSorensen.table(tab, d0 = d0,
                                 conf.level = conf.level, check.table = FALSE))
}

#' @describeIn equivTestSorensen S3 method for class "list"
#' @export
equivTestSorensen.list <- function(x, d0 = 1 / (1 + 1.25),
                                   conf.level = 0.95, ...){
  numLists <- length(x)
  lstNams <- names(x)
  equivTests <- lapply(2:numLists, function(iLst1) {
    oneVsOthers <- lapply(1:(iLst1-1), function(iLst2) {
      return(equivTestSorensen(x[[iLst1]], x[[iLst2]], d0 = d0, conf.level = conf.level,
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
#' @param ontos "character", GO ontologies to analyse. Defaults to \code{c("BP", "CC", "MF")}.
#' @param GOLevels "integer", GO levels to analyse inside each of these GO ontologies.
#' @param ... extra parameters for function \code{crossTabGOIDs4GeneLists} in package \code{equivStandardTest}.
#'
#' @return
#' An object of class "AllEquivSDhtest". It is a list with as many components as GO ontologies have been analysed.
#' Each of these elements is itself a list with as many components as GO levels have been analized.
#' Finally, the elements of these lists are objects as generated by \code{equivTestSorensen.list},
#' i.e., objects of class "equivSDhtestList" containing all pairwise comparisons between the gene lists
#' in argument \code{x}.
#'
#' @examples
#' library(equivStandardTest)
#' # Gene universe:
#' data(humanEntrezIDs)
#' # Gene lists to be explored for enrichment:
#' data(pbtGeneLists)
#' AllEquivTestSorensen(pbtGeneLists,
#'                     geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#'                     ontos = c("MF", "BP"), GOLevels = 4:6)

#'
#' @export
allEquivTestSorensen <- function(x, d0 = 1 / (1 + 1.25), conf.level = 0.95,
                                 ontos = c("BP", "CC", "MF"),
                                 GOLevels = 3:10,
                                 ...)
{
  allOntos <- lapply(ontos, function(onto) {
    thisOnto <- lapply(GOLevels, function(lev) {
      cat("========================================================\n")
      cat("Ontology ", onto, " level ", lev, "\n")
      result <- equivTestSorensen(x, d0 = d0, conf.level = conf.level,
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

