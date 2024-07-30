#' This function builds a cross-tabulation of enriched (TRUE) and non-enriched (FALSE) GO terms vs. gene lists
#'
#' @param x either an object of class "character" (or coerzable to "character") or "list".
#' In the "character" interface, these values should represent Entrez gene (or, in general, feature)
#' identifiers.
#' In the "list" interface, each element of the list must be a "character" vector of
#' Entrez identifiers
#' @param orgPackg A string with the name of the genomic annotation package corresponding to a specific species to be analyzed, which must be previously installed and activated. For more details see \href{../doc/README.html}{README File}.
#' @param geneUniverse character vector containing the universe of genes from where gene lists have been extracted. This vector must be obtained from the annotation package declared in \code{orgPackg}. For more details see \href{../doc/README.html}{README File}.
#' @param onto string describing the ontology. Belongs to c('BP', 'MF', 'CC')
#' @param GOLevel GO level, an integer
# @param restricted Boolean variable to decide how tabulation of GOIDs is performed.
# Unrestricted tabulation crosses _all_ GO Terms located at the level indicated by `GOLevel` with the two GOIDs lists
# Restricted tabulation crosses only terms from the selected GO level that are _common to ancestor terms of either list_.
# That is, if one term in the selected GO level is not an ancestor of at least one of the gene list most specific GO terms
# it is excluded from the GO Level's terms because it is impossible that it appears as being enriched.
#' @param pAdjustMeth string describing the adjust method. Belongs to c('BH', 'BY', 'Bonf')
#' @param pvalCutoff adjusted pvalue cutoff on enrichment tests to report
#' @param qvalCutoff qvalue cutoff on enrichment tests to report as significant.
#' Tests must pass i) pvalueCutoff on unadjusted pvalues, ii) pvalueCutoff on adjusted pvalues and
#' iii) qvalueCutoff on qvalues to be reported
#' @param parallel Logical. Only in "list" interface. Defaults to FALSE but put it at TRUE for parallel computation
#' @param nOfCores Number of cores for parallel computations. Only in "list" interface
#' @param ... Additional parameters
#' 
#' @return In the "character" interface, a length k vector of TRUE/FALSE values
#' corresponding to enrichment or not, where k stands for the total number of GO terms at level
#' 'GOLev' in ontology 'onto'.
#' In the "list" interface, a boolean matrix of TRUE/FALSE values indicating enrichment or not,
#' with k rows and s columns, where k corresponds to the total number of GO terms at level 'GOLev'
#' in ontology 'onto' and s corresponds to the length of "list" 'x'.
#'

#' @details
#' The arguments 'parallel' and 'nOfCores' are ignored in the 'default' and "character" interfaces
#' because (in the present implementation) parallelisation is only applied to repeated calls
#' to function 'clusterProfiler::enrichGO' which, in turn, does not provide for the possibility
#' of parallelisation. They only apply to the "list" interface.
#' 

#' @importFrom org.Hs.eg.db org.Hs.eg.db 
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
#' # Computing the cross table:
#' enrichd <- enrichedIn(allOncoGeneLists[["Vogelstein"]],
#'                       geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#'                       onto = "MF", GOLevel = 6)
#' enrichd
#'
#' # Cross table of enriched GO terms (GO ontology MF, level 6) for all gene 
#' # lists in 'allOncoGeneLists':
#' enrichedAllOncoMF.6 <- enrichedIn(allOncoGeneLists,
#'                           geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#'                           onto = "MF", GOLevel = 6)
#' enrichedAllOncoMF.6

#' @export
enrichedIn <- function(x, ...) {
  UseMethod("enrichedIn")
}

#' @describeIn enrichedIn S3 default method
#' @export
enrichedIn.default <- function (x, geneUniverse, orgPackg,
                                onto, GOLevel,
                                pAdjustMeth = "BH", pvalCutoff = 0.01, qvalCutoff = 0.05,
                                parallel = FALSE,
                                nOfCores = 1, ...){
  if (!requireNamespace(orgPackg, quietly = TRUE)) {
    stop(paste("Genomic annotation of the organism to analyse is in package", orgPackg, ". Please, install this package before to use this function."),
         call. = FALSE)
  }
  allGOIDs <- goProfiles::getGOLevel(onto, GOLevel)
  enriched <- clusterProfiler::enrichGO(gene = as.character(x),
                                        universe = geneUniverse, OrgDb = orgPackg,
                                        ont = onto, pAdjustMethod = pAdjustMeth,
                                        pvalueCutoff = pvalCutoff, qvalueCutoff = qvalCutoff)
  GOIDs <- as.character(as.data.frame(enriched)$ID)
  return(is.element(allGOIDs, GOIDs))
}

#' @describeIn enrichedIn S3 method for class "character"
#' @export
enrichedIn.character <- function (x, geneUniverse, orgPackg,
                                  onto, GOLevel, #restricted = FALSE,
                                  pAdjustMeth = "BH", pvalCutoff = 0.01, qvalCutoff = 0.05,
                                  parallel = FALSE,
                                  nOfCores = 1, ...){
  if (!requireNamespace(orgPackg, quietly = TRUE)) {
    stop(paste("Genomic annotation of the organism to analyse is in package", orgPackg, ". Please, install this package before to use this function."),
         call. = FALSE)
  }
  allGOIDs <- goProfiles::getGOLevel(onto, GOLevel)
  enriched <- clusterProfiler::enrichGO(gene = x,
                                        universe = geneUniverse, OrgDb = orgPackg,
                                        ont = onto, pAdjustMethod = pAdjustMeth,
                                        pvalueCutoff = pvalCutoff, qvalueCutoff = qvalCutoff)
  GOIDs <- as.character(as.data.frame(enriched)$ID)
  result <- is.element(allGOIDs, GOIDs)
  names(result) <- allGOIDs
  return(result)
}

#' @describeIn enrichedIn S3 method for class "list"
#' @importFrom parallel detectCores makeCluster clusterExport clusterEvalQ parSapply stopCluster
#' @export
enrichedIn.list <- function (x, geneUniverse, orgPackg,
                             onto, GOLevel, #restricted = FALSE,
                             pAdjustMeth = "BH", pvalCutoff = 0.01, qvalCutoff = 0.05,
                             parallel = FALSE,
                             nOfCores = min(detectCores() - 1, length(x)), ...){
  if (!requireNamespace(orgPackg, quietly = TRUE)) {
    stop(paste("Genomic annotation of the organism to analyse is in package", orgPackg, ". Please, install this package before to use this function."),
         call. = FALSE)
  }
  allGOIDs <- goProfiles::getGOLevel(onto, GOLevel)
  lenAllGOIDs <- length(allGOIDs)
  lenGeneLists <- length(x)
  if (parallel) {
    on.exit(stopCluster(cl))
    if (.Platform$OS.type == "windows") {
      cl <- makeCluster(nOfCores)
      clusterExport(cl, c("geneUniverse", "orgPackg",
                          "onto", "GOLevel", #"restricted",
                          "pAdjustMeth", "pvalCutoff", "qvalCutoff"),
                    envir = environment())
      clusterEvalQ(cl, {
        library(orgPackg, character.only = TRUE)
      })
    } else {
      cl <- makeCluster(nOfCores, type = "FORK")
    }
    result <- parSapply(cl, seq_len(lenGeneLists), function(iList, ...) {
      enriched <- clusterProfiler::enrichGO(gene = x[[iList]],
                                            universe = geneUniverse, OrgDb = orgPackg,
                                            ont = onto, pAdjustMethod = pAdjustMeth,
                                            pvalueCutoff = pvalCutoff, qvalueCutoff = qvalCutoff)
      GOIDs <- as.character(as.data.frame(enriched)$ID)
      return(is.element(allGOIDs, GOIDs))
    })
  } else {
    result <- vapply(seq_len(lenGeneLists), function(iList, ...) {
      enriched <- clusterProfiler::enrichGO(gene = x[[iList]],
                                            universe = geneUniverse, OrgDb = orgPackg,
                                            ont = onto, pAdjustMethod = pAdjustMeth,
                                            pvalueCutoff = pvalCutoff, qvalueCutoff = qvalCutoff)
      GOIDs <- as.character(as.data.frame(enriched)$ID)
      return(is.element(allGOIDs, GOIDs))
    }, FUN.VALUE = logical(lenAllGOIDs))
  }
  rownames(result) <- allGOIDs
  colnames(result) <- names(x)
  return(result)
}
