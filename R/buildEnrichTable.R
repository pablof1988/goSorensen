
#' Creates a 2x2 enrichment contingency table from two gene lists, or all pairwise contingency
#' tables for a "list" of gene lists.
#'
#' @param x either an object of class "character" (or coerzable to "character") representing a
#' vector of gene identifiers or an object of class "list". In this second case, each element
#' of the list must be a "character" vector of gene identifiers. Then, all pairwise
#' contingency tables between these gene lists are built.
#' @param y an object of class "character" (or coerzable to "character") representing a vector of gene identifiers.
#' @param listNames a character(2) with the gene lists names originating the cross-tabulated enrichment frequencies.
#' @param check.table Logical The resulting table must be checked. Defaults to TRUE.
#' @param geneUniverse character vector containing all genes from where geneLists have been extracted.
#' @param orgPackg A string with the name of the annotation package.
#' @param onto string describing the ontology. Either "BP", "MF" or "CC".
#' @param GOLevel An integer, the GO ontology level.
#' @param restricted Logical variable to decide how tabulation of GOIDs is performed. Defaults to FALSE.
#' See the details section.
#' @param pAdjustMeth string describing the adjust method, either "BH", "BY" or "Bonf", defaults to 'BH'.
#' @param pvalCutoff A numeric value. Defaults to 0.01.
#' @param qvalCutoff A numeric value. Defaults to 0.05.
#' @param parallel Logical. Defaults to FALSE but put it at TRUE for parallel computation.
#' @param nOfCores Number of cores for parallel computations. Only in "list" interface.
#' @param ... Additional parameters for internal use (not used for the moment)
#'
#' @return in the "character" interface, an object of class "table" is returned.
#' It represents a 2x2 contingency table interpretable as the cross-tabulation of
#' the enriched GO items in two gene lists: "Number of enriched items in list 1
#' (TRUE, FALSE)" x "Number of enriched items in list 2 (TRUE, FALSE)".
#' In the "list" interface, the result is an object of class "tableList" with all
#' pairwise tables.
#' Class "tableList" corresponds to objects representing all mutual enrichment contingency tables
#' generated in a pairwise fashion:
#' Given gene lists (i.e. "character" vectors of gene identifiers) l1, l2, ..., lk,
#' an object of class "tableList" is a list of lists of
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
#' $lk$l1$t(k,1), $lk$l2$t(k,2), ..., $lk$l(k-1)t(K,k-1)
#'
#' @details
#' Unrestricted tabulation crosses _all_ GO Terms located at the level indicated by `GOLev`
#' with the two GOIDs lists. Restricted tabulation crosses only terms from the selected GO level
#' that are _common to ancestor terms of either list_.
#' That is, if one term in the selected GO level is not an ancestor of at least one of the gene list
#' most specific GO terms it is excluded from the GO Level's terms because it is impossible that it
#' appears as being enriched.

#' @examples
#' # Gene universe:
#' data(humanEntrezIDs)
#' # Gene lists to be explored for enrichment:
#' data(allOncoGeneLists)
#' ?allOncoGeneLists
#' # Table of mutual GO node enrichment between gene lists Vogelstein and sanger,
#' # for ontology MF at GO level 6 (only first 50 genes, to improve speed).
#' vog.VS.sang <- buildEnrichTable(allOncoGeneLists[["Vogelstein"]][seq_len(50)],
#'                                 allOncoGeneLists[["sanger"]][seq_len(50)],
#'                                 geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#'                                 onto = "MF", GOLevel = 6, listNames = c("Vogelstein", "sanger"))
#' vog.VS.sang
#' # This is an inadequate table for Sorensen-Dice computations:
#' equivTestSorensen(vog.VS.sang)
#' # This sometimes happens, due too small gene lists or due to poor incidence
#' # of enrichment.
#' #
#' # In fact, the complete gene lists generate a much interesting contingency table:
#' # vog.VS.sang <- buildEnrichTable(allOncoGeneLists[["Vogelstein"]],
#' #                                 allOncoGeneLists[["sanger"]],
#' #                                 geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#' #                                 onto = "MF", GOLevel = 6, listNames = c("Vogelstein", "sanger"))
#' # vog.VS.sang
#' # equivTestSorensen(vog.VS.sang)

#' @export
buildEnrichTable <- function(x, ...) {
  UseMethod("buildEnrichTable")
}

#' @describeIn buildEnrichTable S3 default method
#' @export
buildEnrichTable.default <- function(x, y,
                                     listNames = c("gene.list1", "gene.list2"),
                                     check.table = TRUE, geneUniverse, orgPackg, onto, GOLevel,
                                     restricted = FALSE,
                                     pAdjustMeth = "BH", pvalCutoff = 0.01, qvalCutoff = 0.05, ...)
{
  buildEnrichTable.character(as.character(x), as.character(y), listNames, check.table,
                             geneUniverse, orgPackg, onto, GOLevel,
                             restricted,
                             pAdjustMeth, pvalCutoff, qvalCutoff, ...)
}

#' @describeIn buildEnrichTable S3 method for class "character"
#' @export
buildEnrichTable.character <- function(x, y, listNames = c("gene.list1", "gene.list2"),
                                       check.table = TRUE, geneUniverse, orgPackg, onto, GOLevel,
                                       restricted = FALSE,
                                       pAdjustMeth = "BH", pvalCutoff = 0.01, qvalCutoff = 0.05, ...) {
  stopifnot(
    "Argument 'y' is missing, 'x' and 'y' must be 'character' vectors of valid gene identifiers" =
      !missing(y))
  stopifnot("Arguments 'x' and 'y' must be 'character' vectors of valid gene identifiers" =
              is.character(y))
  stopifnot("Argument 'check.table' must be logical" = is.logical(check.table))
  stopifnot("Arguments 'pValCutoff' and 'qValCutoff' must be numeric between 0 and 1" =
              (0 < pvalCutoff) && (pvalCutoff < 1) && (0 < qvalCutoff) && (qvalCutoff < 1))
  tab <- crossTabGOIDs4GeneLists (genelist1 = x, genelist2 = y,
                                  geneUniverse, orgPackg, onto, GOLevel,
                                  restricted,
                                  pAdjustMeth, pvalCutoff, qvalCutoff, ...)
  if (!all(dim(tab) == 2)) {
    tab <- completeTable(tab)
  }
  tab <- tab[c(2,1),c(2,1)]
  if (check.table){
    nice2x2Table.table(tab)
  }
  dimnames(tab) <- list(c(TRUE, FALSE), c(TRUE, FALSE))
  names(dimnames(tab)) <- paste0("Enriched in ", listNames)
  attr(tab, "onto") <- onto
  attr(tab, "GOLevel") <- GOLevel
  return(tab)
}

#' @describeIn buildEnrichTable S3 method for class "list"
#' @importFrom parallel detectCores makeCluster clusterExport clusterEvalQ parLapply stopCluster
#' @export
buildEnrichTable.list <- function(x,
                                  check.table = TRUE, geneUniverse, orgPackg, onto, GOLevel,
                                  restricted = FALSE,
                                  pAdjustMeth = "BH", pvalCutoff = 0.01, qvalCutoff = 0.05,
                                  parallel = FALSE,
                                  nOfCores = min(detectCores() - 1, length(x) - 1),
                                  ...)
{
  lstNams <- names(x)
  numLists <- length(x)
  if (parallel) {
    on.exit(stopCluster(cl))
    if (.Platform$OS.type == "windows") {
      cl <- makeCluster(nOfCores)
      clusterExport(cl, c("lstNams", "geneUniverse", "orgPackg",
                          "onto", "GOLevel", "restricted",
                          "pAdjustMeth", "pvalCutoff", "qvalCutoff"),
                    envir = environment())
      clusterEvalQ(cl, {
        library(clusterProfiler)
        library(goProfiles)
        library(GO.db)
        library(orgPackg, character.only = TRUE)
      })
    } else {
      cl <- makeCluster(nOfCores, type = "FORK")
    }

    allTables <- parLapply(cl, seq.int(2,numLists), function(iLst1, ...) {
      oneVsOthers <- lapply(seq_len(iLst1-1), function(iLst2, ...) {
        return(buildEnrichTable.character(x[[iLst1]], x[[iLst2]],
                                          listNames = lstNams[c(iLst1, iLst2)],
                                          check.table = check.table,
                                          geneUniverse = geneUniverse, orgPackg = orgPackg,
                                          onto = onto, GOLevel = GOLevel,
                                          restricted = restricted,
                                          pAdjustMeth = pAdjustMeth,
                                          pvalCutoff = pvalCutoff, qvalCutoff = qvalCutoff, ...))
      })
      names(oneVsOthers) <- lstNams[seq_len(iLst1-1)]
      return(oneVsOthers)
    })
  } else {
    allTables <- lapply(seq.int(2,numLists), function(iLst1, ...) {
      oneVsOthers <- lapply(seq_len(iLst1-1), function(iLst2, ...) {
        return(buildEnrichTable.character(x[[iLst1]], x[[iLst2]],
                                          listNames = lstNams[c(iLst1, iLst2)],
                                          check.table = check.table,
                                          geneUniverse = geneUniverse, orgPackg = orgPackg,
                                          onto = onto, GOLevel = GOLevel,
                                          restricted = restricted,
                                          pAdjustMeth = pAdjustMeth,
                                          pvalCutoff = pvalCutoff, qvalCutoff = qvalCutoff, ...))
      })
      names(oneVsOthers) <- lstNams[seq_len(iLst1-1)]
      return(oneVsOthers)
    })
  }
  attr(allTables, "onto") <- onto
  attr(allTables, "GOLevel") <- GOLevel
  names(allTables) <- lstNams[seq.int(2,numLists)]
  class(allTables) <- c("tableList", "list")
  return(allTables)
}

#' Iterate \code{buildEnrichTable} along the specified GO ontologies and GO levels
#'
#' @param x object of class "list". Each of its elements must be a "character" vector of gene
#' identifiers. Then all pairwise contingency tables of joint enrichment are built between
#' these gene lists, iterating the process for all specified GO ontologies and GO levels.
#' @param check.table Boolean. If TRUE (default), all resulting tables are checked by means
#' of function \code{nice2x2Table}.
#' @param ontos "character", GO ontologies to analyse. Defaults to \code{c("BP", "CC", "MF")}.
#' @param GOLevels "integer", GO levels to analyse inside each of these GO ontologies.
#' @param trace Logical. If TRUE (default), the (usually very time consuming) process of function
#' \code{allBuildEnrichTable} is traced along the specified GO ontologies and levels.
#' @param ... extra parameters for function \code{buildEnrichTable}.
#'
#' @return
#' An object of class "allTableList". It is a list with as many components as GO ontologies have been
#' analysed.
#' Each of these elements is itself a list with as many components as GO levels have been analised.
#' Finally, the elements of these lists are objects as generated by \code{buildEnrichTable.list},
#' i.e., objects of class "tableList" containing all pairwise contingency tables of mutual enrichmen
#' between the gene lists in argument \code{x}.
#'
#' @examples
#' # This example is extremely time consuming, it scans two GO ontologies and three
#' # GO levels inside them to obtain the contingency tables of mutual enrichment.
#' # Gene universe:
#' # data(humanEntrezIDs)
#' # Gene lists to be explored for enrichment:
#' # data(pbtGeneLists)
#' # allBuildEnrichTable(pbtGeneLists,
#' #                     geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#' #                     ontos = c("MF", "BP"), GOLevels = seq.int(4,6))
#'
#' @export
allBuildEnrichTable <- function(x, check.table = TRUE,
                                ontos = c("BP", "CC", "MF"),
                                GOLevels = seq.int(3,10),
                                trace = TRUE,
                                ...)
{
  allOntos <- lapply(ontos, function(onto) {
    if (trace) {
      cat("\nBuilding contingency tables for ontology ", onto, "\n")
    }
    thisOnto <- lapply(GOLevels, function(lev) {
      if (trace) {
        cat("\n Level ", lev, "\n")
      }
      result <- buildEnrichTable(x, check.table = check.table,
                                 onto = onto, GOLevel = lev, ...)
      return(result)
    })
    names(thisOnto) <- paste("level", GOLevels)
    return(thisOnto)
  })
  names(allOntos) <- ontos
  class(allOntos) <- c("allTableList", "list")
  return(allOntos)
}

