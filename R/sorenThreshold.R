#' For a given level (2, 3, ...) in a GO ontology (BP, MF or CC), compute the equivalence threshold
#'   dissimilarity matrix.
#'
#' @param x either an object of class "list" or class "tableList".
#' See the details section for more information.
#' @param GOLevel integer (2, 3, ...) level of a GO ontology where the GO profiles are built
#' @param onto character, GO ontology ("BP", "MF" or "CC") under consideration
#' @param orgPackg A string with the name of the genomic annotation package corresponding to a specific species to be analyzed, which must be previously installed and activated. For more details see \href{../doc/README.html}{README File}.
#' @param geneUniverse character vector containing the universe of genes from where geneLists have been extracted. This vector must be extracted from the annotation package declared in \code{orgPackg}. For more details see \href{../doc/README.html}{README File}.
#' @param boot boolean. If TRUE, the p-values are computed by means
#' of a bootstrap approach instead of the asymptotic normal approach. Defaults to FALSE.
#' @param nboot numeric, number of initially planned bootstrap replicates. Ignored if
#' \code{boot == FALSE}. Defaults to 10000
#' @param boot.seed starting random seed for all bootstrap iterations. Defaults to 6551.
#'   see the details section
#' @param trace boolean, the full process must be traced? Defaults to TRUE
#' @param alpha simultaneous nominal significance level for the equivalence tests to be repeteadly performed,
#'   defaults to 0.05
#' @param precis numerical precision in the iterative search of the equivalence threshold dissimilarities,
#'   defaults to 0.001
#' @param ... additional arguments to \code{buildEnrichTable}
#' @return An object of class "dist", the equivalence threshold dissimilarity matrix based on the
#' Sorensen-Dice dissimilarity.
#' @details
#' If \code{x} is an object of class "list", each of its elements must be a "character" vector of gene
#' identifiers (e.g., ENTREZ). Then all pairwise threshold dissimilarities between these gene lists are obtained.
#'
#' Class "tableList" corresponds to objects representing all mutual enrichment contingency tables
#' generated in a pairwise fashion:
#' Given gene lists l1, l2, ..., lk, an object of class "tableList" (typically constructed by a call to function
#' \code{\link{buildEnrichTable}}) is a list of lists of
#' contingency tables tij generated from each pair of gene lists i and j, with the
#' following structure:
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
#' If \code{x} is an object of class "tableList", the threshold dissimilarity is obtained
#' over each one of these tables.
#'
#' If \code{boot == TRUE}, all series of \code{nboot} bootstrap replicates start from the same random
#' seed, provided by the argument \code{boot.seed}, except if \code{boot == NULL}.
#'
#' Do not confuse the resulting threshold dissimilarity matrix with the Sorensen-Dice dissimilarities
#' computed in each equivalence test.
#'
#' The dimension of the resulting matrix may be less than the number of original gene lists being
#' compared, as the process may not converge for some pairs of gene lists.
#'
#' @importFrom stats as.dist
#' @examples
#' # Gene lists to be explored for enrichment:
#' data(allOncoGeneLists)
#' 
#' # Obtaining ENTREZ identifiers for the gene universe of humans:
#' library(org.Hs.eg.db)
#' humanEntrezIDs <- keys(org.Hs.eg.db, keytype = "ENTREZID")
#' 
#' # This example is quite time consuming:
#' # sorenThreshold(allOncoGeneLists,
#' #                geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
#' # Much faster:
#' # Object \code{cont_all_BP4} of class "tableList" contains all the pairwise contingency
#' # tables of joint enrichment for the gene lists in \code{allOncoGeneLists}, for the BP
#' # GO ontology at level 4:
#' data("cont_all_BP4")
#' sorenThreshold(cont_all_BP4)
#'

#' @export
sorenThreshold <- function(x, ...) {
  UseMethod("sorenThreshold")
}

#' @describeIn sorenThreshold S3 method for class "list"
#' @export
sorenThreshold.list <- function(x, onto, GOLevel, geneUniverse, orgPackg,
                                boot = FALSE, nboot = 10000, boot.seed = 6551,
                                trace = TRUE, alpha = 0.05, precis = 0.001, ...)
{
  if (trace) {
    cat("\n", date(), "  Building all enrichment contingency tables...\n")
  }
  all2x2Tables <- buildEnrichTable(x,
                                   onto = onto, GOLevel = GOLevel,
                                   geneUniverse = geneUniverse, orgPackg = orgPackg,
                                   ...)
  return(sorenThreshold.tableList(all2x2Tables, onto = onto, GOLevel = GOLevel,
                                  boot = boot, nboot = nboot, boot.seed = boot.seed,
                                  trace = trace, alpha = alpha, precis = precis, ...))
}

#' @describeIn sorenThreshold S3 method for class "tableList"
#' @export
sorenThreshold.tableList <- function(x, #onto = NULL, GOLevel = NULL,
                                     #geneUniverse, orgPackg,
                                     boot = FALSE, nboot = 10000, boot.seed = 6551,
                                     trace = TRUE, alpha = 0.05, precis = 0.001, ...)
{
  s <- length(x) + 1
  h <- s * (s - 1) * 0.5
  equivDists <- rep(NA, h)
  dIdxs <- unlist(lapply(seq.int(2,s),
                         function(i) lapply(seq.int(1,(i-1)),
                                            function(j) c(i,j))
  ), recursive = FALSE)
  if (boot) {
    arrTabs <- array(unlist(x), dim = c(2,2,h))
  }
  d <- dSorensen(x)
  d <- d[upper.tri(d)]
  se <- seSorensen(x)
  se <- se[upper.tri(se)]
  finite.se <- se > 0.0
  finite.se[!is.finite(se)] <- FALSE
  h.finite <- sum(finite.se)
  if (h.finite == 0) {
    return(NULL)
  }
  finite.idxs <- integer(h.finite)
  ih.finite <- 0
  for (ih in seq.int(1,h)) {
    if (finite.se[ih]) {
      ih.finite <- ih.finite + 1
      finite.idxs[ih.finite] <- ih
    }
  }
  d <- d[finite.se]
  se <- se[finite.se]
  if (boot) {
    arrTabs <- arrTabs[,,finite.se, drop = FALSE]
  }
  deltaMax <- startingDelta(alpha, d, se, h.finite, boot, nboot, boot.seed, arrTabs)
  delta <- deltaMax

  if (trace) {
    cat("\n", date(), "  Computing all threshold equivalence distances...\n")
  }
  for (ih in seq.int(h.finite,1)) {
    nextStep <- nextEquivDist(d, se, ih, delta, alpha, precis,
                              boot, nboot, boot.seed, arrTabs)
    delta <- nextStep$delta
    iDelta <- nextStep$iDelta
    equivDists[finite.idxs[iDelta]] <- delta
    se[iDelta] <- NA
  }
  distMat <- matrix(NA, nrow = s, ncol = s)
  rownames(distMat) <- colnames(distMat) <- c(names(x[[1]]), names(x))
  distMat[upper.tri(distMat)] <- equivDists
  distMat[lower.tri(distMat)] <- t(distMat)[lower.tri(distMat)]

  validToPlot <- apply(distMat, 1, function(thisRow) !all(is.na(thisRow)))
  distMat <- distMat[validToPlot, validToPlot]
  distMat[is.na(distMat)] <- 1.0
  validToPlot <- apply(distMat, 1, function(thisRow) !all(thisRow >= 1.0))
  distMat <- distMat[validToPlot, validToPlot]

  if (nrow(distMat) <= 1) {
    if (trace) {
      cat("\n", date(), "Invalid threshold distance matrix\n")
    }
    return(NULL)
  }
  distMat[distMat > 1.0] <- 1.0
  diag(distMat) <- 0.0
  distMat <- as.dist(t(distMat))
  attr(distMat, "dist.method") <- "equivalence_threshold"
  attr(distMat, "all2x2Tables") <- x
  attr(distMat, "onto") <- attr(x, "onto")
  attr(distMat, "GOLevel") <- attr(x, "GOLevel")
  return(distMat)
}

startingDelta <- function(alpha, d, se, h, boot, nboot, boot.seed, arrTabs) {
  alphas.holm <- alpha / seq.int(h,1)
  delta <- max(c(1, qnorm(1 - alpha / h)) %*% t(matrix(c(d, se), ncol = 2)))
  incDelta <- delta * 0.1
  repeat {
    p.vals <- pvals(delta, d, se, boot, nboot, boot.seed, arrTabs)
    p.order <- order(p.vals)
    p.sort <- p.vals[p.order][seq.int(1,h)]
    if (all(p.sort <= alphas.holm)) {
      return(delta)
    } else {
      delta <- delta + incDelta
    }
  }
}

nextEquivDist <- function(d, se, h, delta, alpha, precis, boot,
                          nboot, boot.seed, arrTabs) {
  incDelta <- 0.1 * delta
  alphas.holm <- alpha / seq.int(h,1)
  repeat {
    p.vals <- pvals(delta, d, se, boot, nboot, boot.seed, arrTabs)
    p.order <- order(p.vals)
    p.sort <- p.vals[p.order][seq.int(1,h)]
    if (all(p.sort <= alphas.holm)) {
      if (incDelta < precis) {
        return(list(delta = delta, iDelta = p.order[h]))
      } else {
        incDelta <- incDelta * 0.5
        delta <- delta - incDelta
      }
    } else {
      delta <- delta + incDelta
      incDelta <- 0.5 * incDelta
      delta <- delta - incDelta
    }
  }
}

pvals <- function(delta, d, se, boot, nboot, boot.seed, arrTabs){
  if (boot) {
    if (!is.null(boot.seed)) {
      set.seed(boot.seed)
    }
    p.vals <- unlist(vapply(seq.int(1,length(se)), function(i) {
      if (is.na(se[i])) {
        return(NA)
      } else {
        di <- d[i]
        return(pboot((di - delta) / se[i], nboot, arrTabs[,,i], di))
      }
    }, FUN.VALUE = numeric(1)))
  } else {
    p.vals <- ifelse(is.na(se), NA, pnorm((d - delta) / se))
  }
}

pboot <- function(x, nboot, tab, d) {
  n <- sum(tab)
  nu <- sum(tab[seq.int(1,3)])
  pTab <- tab / n
  bootTabs <- rmultinom(nboot, size = n, prob = pTab)
  tStats <- apply(bootTabs, 2, tStatSorensen, dis = d)
  tStats <- tStats[is.finite(tStats)]
  len.tStats <- length(tStats)
  return((sum(tStats <= x) + 1) / (len.tStats + 1))
}

tStatSorensen <- function(xBoot, dis) {
  nu <- sum(xBoot[seq_len(3)])
  dBoot <- (xBoot[2] + xBoot[3]) / (nu + xBoot[1])
  p11 <- xBoot[1] / nu
  p11plus <- 1 + p11
  se <- 2 * sqrt(p11 * (1 - p11) / (nu - 1)) / (p11plus * p11plus)
  return((dBoot - dis) / se)
}

