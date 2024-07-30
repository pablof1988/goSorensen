
#' Creates a 2x2 enrichment contingency table from two gene lists, or all pairwise contingency
#' tables for a "list" of gene lists.
#'
#' @param x either an object of class "character" (or coerzable to "character") representing a
#' vector of gene identifiers (e.g., ENTREZ) or an object of class "list". In this second case, each
#' element of the list must be a "character" vector of gene identifiers (e.g., ENTREZ). Then, all pairwise
#' contingency tables between these gene lists are built.
#' @param y an object of class "character" (or coerzable to "character") representing a vector
#' of gene identifiers (e.g., ENTREZ).
#' @param listNames a character(2) with the gene lists names originating the cross-tabulated
#' enrichment frequencies. Only in the "character" or default interface.
#' @param check.table Logical The resulting table must be checked. Defaults to TRUE.
#' @param orgPackg A string with the name of the genomic annotation package corresponding to a specific species to be analyzed, which must be previously installed and activated. For more details see \href{../doc/README.html}{README File}.
#' @param geneUniverse character vector containing the universe of genes from where gene lists have been extracted. This vector must be obtained from the annotation package declared in \code{orgPackg}. For more details see \href{../doc/README.html}{README File}.
#' @param onto string describing the ontology. Either "BP", "MF" or "CC".
#' @param GOLevel An integer, the GO ontology level.
#' @param showEnrichedIn Boolean. If TRUE (default), the cross-table of enriched and non-enriched GO terms vs  Gene Lists names (obtained from the \code{enrichedIn} function) is automatically saved in the Global Environment.
#' @param pAdjustMeth string describing the adjust method, either "BH", "BY" or "Bonf", defaults to 'BH'.
#' @param pvalCutoff adjusted pvalue cutoff on enrichment tests to report
#' @param qvalCutoff qvalue cutoff on enrichment tests to report as significant.
#' Tests must pass i) pvalueCutoff on unadjusted pvalues, ii) pvalueCutoff on adjusted pvalues and
#' iii) qvalueCutoff on qvalues to be reported
#' @param parallel Logical. Defaults to FALSE but put it at TRUE for parallel computation.
#' @param nOfCores Number of cores for parallel computations. Only in "list" interface.
#' @param ... Additional parameters for internal use (not used for the moment)
#'
#' @return in the "character" interface, an object of class "table".
#' It represents a 2x2 contingency table, the cross-tabulation of
#' the enriched GO terms in two gene lists: "Number of enriched GO terms in list 1
#' (TRUE, FALSE)" x "Number of enriched Go terms in list 2 (TRUE, FALSE)".
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
#' The arguments 'parallel' and 'nOfCores' are ignored in the 'default' and "character" interfaces,
#' but included for possible future developments; they only apply to the "list" interface.
#' In the "list" interface, 'parallel' defaults to FALSE but there is the possibility of some time
#' saving when the number of gene lists (the length of 'x' in the "list" interface) is high.
#' The trade off between the time spent initializing parallel computing and the possible time gain
#' due to parallelization must be considered in each application and computer.
#' 

#' @examples
#' # Obtaining ENTREZ identifiers for the gene universe of humans:
#' library(org.Hs.eg.db)
#' humanEntrezIDs <- keys(org.Hs.eg.db, keytype = "ENTREZID")
#' 
#' # Gene lists to be explored for enrichment:
#' data(allOncoGeneLists)
#' ?allOncoGeneLists
#' 
#' # Table of joint GO term enrichment between gene lists Vogelstein and sanger,
#' # for ontology MF at GO level 6.
#' vog.VS.sang <- buildEnrichTable(allOncoGeneLists[["Vogelstein"]],
#'                                 allOncoGeneLists[["sanger"]],
#'                                 geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#'                                 onto = "MF", GOLevel = 6, listNames = c("Vogelstein", "sanger"))
#' vog.VS.sang
#' # All tables of mutual enrichment:
#' all.tabs <- buildEnrichTable(allOncoGeneLists,
#'                              geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#'                              onto = "MF", GOLevel = 6)
#' all.tabs$waldman

#' @export
buildEnrichTable <- function(x, ...) {
  UseMethod("buildEnrichTable")
}

#' @describeIn buildEnrichTable S3 default method
#' @export
buildEnrichTable.default <- function(x, y,
                                     listNames = c("gene.list1", "gene.list2"),
                                     check.table = TRUE, geneUniverse, orgPackg, onto, GOLevel, showEnrichedIn = TRUE,
                                     pAdjustMeth = "BH", pvalCutoff = 0.01, qvalCutoff = 0.05,
                                     parallel = FALSE, nOfCores = 1, ...){
  if (!requireNamespace(orgPackg, quietly = TRUE)) {
    stop(paste("Genomic annotation of the organism to analyse is in package", orgPackg, ". Please, install this package before to use this function."),
         call. = FALSE)
  }
  buildEnrichTable.character(as.character(x), as.character(y), listNames, check.table,
                             geneUniverse, orgPackg, onto, GOLevel,
                             pAdjustMeth, pvalCutoff, qvalCutoff,
                             parallel, nOfCores, ...)
}

#' @describeIn buildEnrichTable S3 method for class "character"
#' @export
buildEnrichTable.character <- function(x, y, listNames = c("gene.list1", "gene.list2"),
                                      geneUniverse, orgPackg, onto, GOLevel, showEnrichedIn = TRUE, 
                                      check.table = TRUE, pAdjustMeth = "BH", pvalCutoff = 0.01, 
                                      qvalCutoff = 0.05, parallel = FALSE, nOfCores = 1, ...) {
  if (!requireNamespace(orgPackg, quietly = TRUE)) {
    stop(paste("Genomic annotation of the organism to analyse is in package", orgPackg, ". Please, install this package before to use this function."),
         call. = FALSE)
  }
  stopifnot(
    "Argument 'y' is missing, 'x' and 'y' must be 'character' vectors of valid gene identifiers" =
      !missing(y))
  stopifnot("Arguments 'x' and 'y' must be 'character' vectors of valid gene identifiers" =
              is.character(y))
  stopifnot("Argument 'check.table' must be logical" = is.logical(check.table))
  stopifnot("Arguments 'pValCutoff' and 'qValCutoff' must be numeric between 0 and 1" =
              (0 < pvalCutoff) && (pvalCutoff < 1) && (0 < qvalCutoff) && (qvalCutoff < 1))
  enrich1 <- enrichedIn(x, geneUniverse, orgPackg, onto, GOLevel,
                        pAdjustMeth, pvalCutoff, qvalCutoff)
  enrich2 <- enrichedIn(y, geneUniverse, orgPackg, onto, GOLevel,
                        pAdjustMeth, pvalCutoff, qvalCutoff)
  
  if (is.null(listNames)) {
    dnnNames <- paste0("Enriched_in_", c("list_1", "list_2"))
  }else{
    dnnNames <- paste0("Enriched_in_", listNames)
  }
  tab <- table(enrich1, enrich2, dnn=dnnNames)
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
  if(showEnrichedIn){
    enrichCross <- cbind(enrich1, enrich2)
    colnames(enrichCross) <- listNames
    assign(paste0("enrichedIn_", onto, GOLevel), enrichCross, envir = .GlobalEnv)
  }
  return(tab)
}

#' @describeIn buildEnrichTable S3 method for class "list"
# @importFrom parallel detectCores makeCluster clusterExport clusterEvalQ parLapply stopCluster
#' @export
buildEnrichTable.list <- function(x,
                                 check.table = TRUE, geneUniverse, orgPackg, onto, GOLevel,
                                 showEnrichedIn = TRUE, pAdjustMeth = "BH", pvalCutoff = 0.01, 
                                 qvalCutoff = 0.05, parallel = FALSE,
                                 nOfCores = min(detectCores() - 1, length(x) - 1),
                                 ...){
  if (!requireNamespace(orgPackg, quietly = TRUE)) {
    stop(paste("Genomic annotation of the organism to analyse is in package", orgPackg, ". Please, install this package before to use this function."),
         call. = FALSE)
  }
  lstNams <- names(x)
  numLists <- length(x)
  allEnrichs <- enrichedIn(x, geneUniverse, orgPackg, onto, GOLevel, #restricted = FALSE,
                           pAdjustMeth, pvalCutoff, qvalCutoff, parallel, nOfCores)
  allTables <- lapply(seq.int(2,numLists), function(iLst1, ...) {
    oneVsOthers <- lapply(seq_len(iLst1-1), function(iLst2, ...) {
      tab <- table(allEnrichs[,iLst1], allEnrichs[,iLst2],
                   dnn = paste0("Enriched_in_", lstNams[c(iLst1, iLst2)]))
      if (!all(dim(tab) == 2)) {
        tab <- completeTable(tab)
      }
      tab <- tab[c(2,1),c(2,1)]
      if (check.table){
        nice2x2Table.table(tab)
      }
      dimnames(tab) <- list(c(TRUE, FALSE), c(TRUE, FALSE))
      names(dimnames(tab)) <- paste0("Enriched in ", lstNams[c(iLst1, iLst2)])
      attr(tab, "onto") <- onto
      attr(tab, "GOLevel") <- GOLevel
      return(tab)
    })
    names(oneVsOthers) <- lstNams[seq_len(iLst1-1)]
    return(oneVsOthers)
  })

  attr(allTables, "onto") <- onto
  attr(allTables, "GOLevel") <- GOLevel
  names(allTables) <- lstNams[seq.int(2,numLists)]
  class(allTables) <- c("tableList", "list")
  if(showEnrichedIn){
    assign(paste0("enrichedIn_", onto, GOLevel), allEnrichs, envir = .GlobalEnv)
  }
  return(allTables)
}

