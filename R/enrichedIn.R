#' This function builds a cross-tabulation of enriched (TRUE) and non-enriched
#' (FALSE) GO terms vs. gene lists
#'
#' @param x Either an object of class \code{"character"} (or coercible to
#'   \code{"character"}) or \code{"list"}.
#'   In the \code{"character"} interface, these values should represent Entrez
#'   gene identifiers (or, more generally, feature identifiers).
#'   In the \code{"list"} interface, each element of the list must be a
#'   character vector of identifiers.
#' @param ... Additional parameters passed to internal method implementations.
#' @export
enrichedIn <- function(x, ...) {
  UseMethod("enrichedIn")
}
#' @param orgPackg A string with the name of the genomic annotation package
#' corresponding to a specific species to be analyzed, which must be previously
#' installed and activated.
#' For more details, refer to vignette
#' \href{../doc/goSorensen_Introduction.html}{goSorensen_Introduction}.
#' @param geneUniverse Character vector containing the universe of genes from
#' which the gene lists were extracted.
#'   This vector must be obtained from the annotation package declared in
#'   \code{orgPackg}.
#'   For more details, refer to vignette
#'   \href{../doc/goSorensen_Introduction.html}{goSorensen_Introduction}.
#' @param onto string describing the ontology. Belongs to c('BP', 'MF', 'CC')
#' @param GOLevel Integer GO level to analyze within the selected ontology.
#'   If \code{NULL}, the analysis is performed without restricting the
#'   enrichment computation to a specific GO level.
#' @param pAdjustMeth string describing the adjust method. Belongs to c('BH',
#' 'BY', 'Bonf')
#' @param pvalCutoff adjusted pvalue cutoff on enrichment tests to report
#' @param qvalCutoff qvalue cutoff on enrichment tests to report as significant.
#' Tests must pass i) pvalueCutoff on unadjusted pvalues, ii) pvalueCutoff on
#' adjusted pvalues and iii) qvalueCutoff on qvalues to be reported
#' @param parallel Logical. Only in "list" interface. Defaults to FALSE, put it
#' at TRUE for parallel computation
#' @param nOfCores Number of cores for parallel computations. Only in "list"
#' interface
#' @param onlyEnriched logical. If TRUE (the default), the returned result only
#' contains those GO terms
#' enriched in almost one of the gene lists
#' @param keyType Character string giving the type of gene identifier used in
#' \code{x}, such as \code{"ENTREZID"} or \code{"SYMBOL"}.
#' @param ... Additional parameters passed to internal method implementations.
#'
#' @return In the "character" interface, a length k vector of TRUE/FALSE values
#' corresponding to enrichment or not of the GO terms at level 'GOLevel' in
#' ontology 'onto'.
#' If 'onlyEnriched' is FALSE, k corresponds to the total number of these GO
#' terms. If 'onlyEnriched'
#' is TRUE (default), k is the number of enriched GO terms (and then all values
#' in the resulting vector are TRUE).
#' In the "list" interface, a logical matrix of TRUE/FALSE values indicating
#' enrichment or not, with k rows and s columns. s is the number of gene lists
#' (the length of list 'x').
#' If 'onlyEnriched' is FALSE, k corresponds to the total number of GO terms at
#' level 'GOLevel' in ontology onto'. If 'onlyEnriched' is TRUE (default), the
#' resulting matrix only contains the k rows corresponding to GO terms enriched
#' in at least one of these s gene lists.
#' In both interfaces ("character" or "list"), the result also has an attribute
#' (\code{nTerms}) with the total number of GO terms at level 'GOLevel' in
#' ontology 'onto'.
#'
#' @details
#' When the function argument \code{onlyEnriched} is FALSE, commonly the result
#' is a sparse but very large object. This function is primarily designed for
#' internal use of function \code{buildEnrichTable}, with argument
#' \code{onlyEnriched} always put at its default TRUE value.
#' Then calls to \code{enrichedIn} result in much more compact objects, in
#' general.
#'
#' Argument \code{parallel} only applies to interface "list". Its default value
#' is "FALSE" and you may consider the trade of between the time spent in
#' initializing parallelization and the possible time gain when parallelizing.
#' It is difficult to establish a general guideline, but parallelizing
#' is only worthwhile when analyzing many gene lists, on the order of 30 or
#' more, although it depends a lot on each processor.
#'
#' AnnotationDbi::select(org.Hs.eg.db, ...)
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
#' # iii) Computing the threshold enrichment matrix directly from gene lists at
#' # GO level 4 and BP ontology:
#' enrichedInBP4 <- enrichedIn(allOncoGeneLists,
#'   geneUniverse = humanEntrezIDs,
#'   orgPackg = "org.Hs.eg.db",
#'   onto = "BP",
#'   GOLevel = 4
#' )
#' enrichedInBP4
#' }
#'
#' # Since running this example may take several minutes, the result has been
#' # pre-computed and is accessible as the following:
#' data(enrichedInBP4)
#' enrichedInBP4
#' # This shortcut applies only to this example; for your own gene-list data,
#' # the computation must be performed explicitly.
#'
#' # For a complete overview of this function's use, see the section 2 of the
#' # vignette 'Introduction to goSorensen'. You can do this by consulting the
#' # general package documentation or by directly running the following code in
#' # the R console:
#' # vignette("goSorensen_Introduction", package = "goSorensen")
#'
#' @describeIn enrichedIn S3 default method
#' @export
enrichedIn.default <- function(x, geneUniverse, orgPackg,
                               onto, GOLevel,
                               pAdjustMeth = "BH", pvalCutoff = 0.01,
                               qvalCutoff = 0.05,
                               parallel = FALSE,
                               nOfCores = 1,
                               onlyEnriched = TRUE,
                               keyType = "ENTREZID", ...) {
  enrichedIn.character(
    as.character(x),
    geneUniverse = geneUniverse,
    orgPackg = orgPackg,
    onto = onto,
    GOLevel = GOLevel,
    pAdjustMeth = pAdjustMeth,
    pvalCutoff = pvalCutoff,
    qvalCutoff = qvalCutoff,
    parallel = parallel,
    nOfCores = nOfCores,
    onlyEnriched = onlyEnriched,
    keyType = keyType,
    ...
  )
}
#' @export
enrichedIn.character <- function(x, geneUniverse, orgPackg,
                                 onto, GOLevel,
                                 pAdjustMeth = "BH", pvalCutoff = 0.01,
                                 qvalCutoff = 0.05,
                                 parallel = FALSE,
                                 nOfCores = 1,
                                 onlyEnriched = TRUE,
                                 keyType = "ENTREZID", ...) {
  if (!requireNamespace(orgPackg, quietly = TRUE)) {
    stop("Genomic annotation of the organism to analyse is in package",
      orgPackg,
      ". Please, install this package before to use this function.",
      call. = FALSE
    )
  }

  orgDbObj <- get(orgPackg, envir = asNamespace(orgPackg))

  enriched <- clusterProfiler::enrichGO(
    gene = x,
    universe = geneUniverse,
    OrgDb = orgPackg,
    ont = onto,
    pAdjustMethod = pAdjustMeth,
    pvalueCutoff = pvalCutoff,
    qvalueCutoff = qvalCutoff,
    keyType = keyType
  )

  goIds <- as.character(as.data.frame(enriched)$ID)

  if (!is.null(GOLevel)) {
    allGoIds <- goProfiles::getGOLevel(onto, GOLevel)
  } else {
    allGoIds <- unique(as.character(enriched@result$ID))
  }

  nGoIds <- length(allGoIds)

  result <- is.element(allGoIds, goIds)
  names(result) <- allGoIds

  if (onlyEnriched) {
    result <- result[result]
  }

  attr(result, "nTerms") <- nGoIds
  return(result)
}
#' @export
enrichedIn.list <- function(x, geneUniverse, orgPackg,
                            onto, GOLevel,
                            pAdjustMeth = "BH", pvalCutoff = 0.01,
                            qvalCutoff = 0.05,
                            parallel = FALSE,
                            nOfCores = min(
                              parallel::detectCores() - 1,
                              length(x)
                            ),
                            onlyEnriched = TRUE,
                            keyType = "ENTREZID", ...) {
  if (!requireNamespace(orgPackg, quietly = TRUE)) {
    stop(
      "Genomic annotation of the organism to analyse is in package",
      orgPackg,
      ". Please, install this package before using this function.",
      call. = FALSE
    )
  }

  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Package 'clusterProfiler' must be installed.", call. = FALSE)
  }

  # if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
  #  stop("Package 'AnnotationDbi' must be installed.", call. = FALSE)
  # }

  orgDbObj <- get(orgPackg, envir = asNamespace(orgPackg))
  lenGeneLists <- length(x)

  if (is.null(names(x))) {
    names(x) <- paste("List", seq_along(x))
  }

  if (!is.null(GOLevel)) {
    allGoIds <- goProfiles::getGOLevel(onto, GOLevel)
    nGoIds <- length(allGoIds)

    if (parallel) {
      on.exit(parallel::stopCluster(cl), add = TRUE)

      if (.Platform$OS.type == "windows") {
        cl <- parallel::makeCluster(nOfCores)
        parallel::clusterExport(
          cl,
          c(
            "x", "geneUniverse", "orgPackg", "onto", "GOLevel",
            "pAdjustMeth", "pvalCutoff", "qvalCutoff", "keyType", "allGoIds"
          ),
          envir = environment()
        )
        parallel::clusterEvalQ(cl, {
          loadNamespace("clusterProfiler")
          NULL
        })
      } else {
        cl <- parallel::makeCluster(nOfCores, type = "FORK")
      }

      result <- parallel::parSapply(
        cl,
        seq_len(lenGeneLists),
        function(iList, ...) {
          enriched <- clusterProfiler::enrichGO(
            gene = x[[iList]],
            universe = geneUniverse,
            OrgDb = orgPackg,
            ont = onto,
            pAdjustMethod = pAdjustMeth,
            pvalueCutoff = pvalCutoff,
            qvalueCutoff = qvalCutoff,
            keyType = keyType,
            ...
          )
          goIds <- as.character(as.data.frame(enriched)$ID)
          is.element(allGoIds, goIds)
        },
        USE.NAMES = FALSE
      )
    } else {
      result <- vapply(
        seq_len(lenGeneLists),
        function(iList, ...) {
          enriched <- clusterProfiler::enrichGO(
            gene = x[[iList]],
            universe = geneUniverse,
            OrgDb = orgPackg,
            ont = onto,
            pAdjustMethod = pAdjustMeth,
            pvalueCutoff = pvalCutoff,
            qvalueCutoff = qvalCutoff,
            keyType = keyType,
            ...
          )
          goIds <- as.character(as.data.frame(enriched)$ID)
          is.element(allGoIds, goIds)
        },
        FUN.VALUE = logical(nGoIds),
        USE.NAMES = FALSE
      )
    }

    rownames(result) <- allGoIds
    colnames(result) <- names(x)

    if (onlyEnriched) {
      result <- result[apply(result, 1, any), , drop = FALSE]
    }

    attr(result, "nTerms") <- nGoIds
    return(result)
  } else {
    getGoData <- function(genes) {
      enriched <- clusterProfiler::enrichGO(
        gene = genes,
        universe = geneUniverse,
        OrgDb = orgDbObj,
        keyType = keyType,
        ont = onto,
        pAdjustMethod = pAdjustMeth,
        pvalueCutoff = pvalCutoff,
        qvalueCutoff = qvalCutoff
      )

      enrichedTable <- enriched@result

      if (is.null(enrichedTable) || nrow(enrichedTable) == 0) {
        return(list(
          enriched = character(0),
          annotated = character(0)
        ))
      }

      goIdsEnriched <- as.character(as.data.frame(enriched)$ID)
      goIdsAnnotated <- unique(as.character(enriched@result$ID))

      list(
        enriched = unique(goIdsEnriched),
        annotated = unique(goIdsAnnotated)
      )
    }

    if (parallel) {
      cl <- parallel::makeCluster(nOfCores)
      on.exit(parallel::stopCluster(cl), add = TRUE)

      parallel::clusterExport(
        cl,
        varlist = c(
          "x", "geneUniverse", "orgDbObj", "onto",
          "pAdjustMeth", "pvalCutoff", "qvalCutoff",
          "keyType", "getGoData"
        ),
        envir = environment()
      )

      parallel::clusterEvalQ(cl, {
        loadNamespace("clusterProfiler")
        NULL
      })

      goDataList <- parallel::parLapply(cl, x, getGoData)
    } else {
      goDataList <- lapply(x, getGoData)
    }

    annotatedList <- lapply(goDataList, function(z) z$annotated)
    enrichedList <- lapply(goDataList, function(z) z$enriched)

    if (length(annotatedList) == 0) {
      result <- matrix(FALSE, nrow = 0, ncol = lenGeneLists)
      colnames(result) <- names(x)
      attr(result, "nTerms") <- 0
      return(result)
    }

    allGoIds <- Reduce(union, annotatedList)
    nGoIds <- length(allGoIds)

    if (nGoIds == 0) {
      result <- matrix(FALSE, nrow = 0, ncol = lenGeneLists)
      colnames(result) <- names(x)
      attr(result, "nTerms") <- 0
      return(result)
    }

    result <- vapply(
      enrichedList,
      function(goIds) is.element(allGoIds, goIds),
      FUN.VALUE = logical(nGoIds)
    )

    rownames(result) <- allGoIds
    colnames(result) <- names(x)

    if (onlyEnriched && nrow(result) > 0) {
      result <- result[apply(result, 1, any), , drop = FALSE]
    }
    attr(result, "nTerms") <- nGoIds
    return(result)
  }
}
