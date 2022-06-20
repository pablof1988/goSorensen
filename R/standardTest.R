#'
#' @importFrom clusterProfiler enrichGO
#' @importFrom goProfiles getGOLevel GOTermsList getAncestorsLst
#' @importFrom stats fisher.test
#' @import GO.db
#' @import org.Hs.eg.db

extractGOIDs <- function (enriched) {
  as.character(as.data.frame(enriched)$ID)
}

#' enrichOnto
#'
#' This function performs standard tests of enrichment from a gene list
#'
#' @param geneList character vector containing a FIRST gene list of entrez IDs
#' @param geneUniverse character vector containing all genes from where geneLists have been extracted
#' @param onto string describing the ontology. Belongs to c('BP', 'MF', 'CC', 'ANY')
#' @param orgPackage A string wih the name of the annotation package
#' @param pAdjustMeth string describing the adjust method. Belongs to c('BH', 'BY', 'Bonf')
#' @param pvalCutoff A numeric value
#' @param qvalCutoff A numeric value

#' @examples
#'  data(kidneyGeneLists)
#'  data(humanEntrezIDs)
#'  gl1 <-kidneyGeneLists[[1]]
#'  anOnto <- 'BP'
#'  enriched <- enrichOnto (geneL=gl1, geneUniverse=humanEntrezIDs, orgPackage='org.Hs.eg.db', onto=anOnto)
#'  GOIDs <- as.character(as.data.frame(enriched)$ID)
enrichOnto <- function (geneList,
                        geneUniverse,
                        orgPackage = 'org.Hs.eg.db',
                        onto = c("BP", "MF", "CC"),
                        pAdjustMeth = "BH",
                        pvalCutoff = 0.01,
                        qvalCutoff = 0.05) {
  enrichGO(gene = geneList, universe = geneUniverse, OrgDb = orgPackage,
           ont = onto,
           pAdjustMethod = pAdjustMeth, pvalueCutoff = pvalCutoff, qvalueCutoff = qvalCutoff,
           readable = TRUE)
}

#' crossTabGOIDsUnrestricted
#'
#' This function performs a crosstabulation between two lists of enriched GOTerms
#' The lists are intended to have been obtained from enrichment analyses performed
#' on two gene lists
#'
#' @param GO1 character vector containing a FIRST list of GO identifiers
#' @param GO2 character vector containing a SECOND gene list GO identifiers
#' @param onto string describing the ontology. Belongs to c('BP', 'MF', 'CC', 'ANY')
#' @param GOLev An integer
#' @param listNames character vector with names of the genelists that generated the
#' enriched GOIDs
#' @examples
#'  data(kidneyEnrichedGOIDs)
#'  GOIDs1 <-kidneyEnrichedGOIDs[[1]]
#'  GOIDs2 <-kidneyEnrichedGOIDs[[2]]
#'  names4lists<- names(kidneyEnrichedGOIDs)[1:2]
#'  anOnto <- 'BP'
#'  GOLev<- 3
#'  crossTabbedGOIds <- crossTabGOIDsUnrestricted (GO1 = GOIDs1, GO2 = GOIDs2, onto = anOnto,
#'  GOLevel =GOLev, listNames=names4lists)
crossTabGOIDsUnrestricted <- function (GO1, GO2, onto, GOLevel, listNames = NULL)
{
  levelIDs <- getGOLevel(onto = onto, level = GOLevel)
  levelInGO1 <- levelIDs %in% GO1
  levelInGO2 <- levelIDs %in% GO2
  if (is.null(listNames)) {
    dnnNames <- paste0("Enriched_in_", c("list_1", "list_2"))
  }else{
    dnnNames <- paste0("Enriched_in_", listNames)
  }
  crossedTable <- table(levelInGO1, levelInGO2, dnn=dnnNames)
}

#' GOIDsInLevel
#'
#' This function extends getGOLevel returning only GO identifiers appearing between GO ancestors of at least one GeneList
#'
#' @param GOLev An integer
#' @param onto string describing the ontology. Belongs to c('BP', 'MF', 'CC', 'ANY')
#' @param geneList1 character vector containing a FIRST gene list of entrez IDs
#' @param geneList2 character vector containing a SECOND gene list of entrez IDs
#' @param orgPackage A string wih the name of the annotation package
#' @param restricted Boolean variable to decide how tabulation is performed.
#' @examples
#'  data(kidneyGeneLists)
#'  GOLev <- 3
#'  anOnto <- 'BP'
#'  gl1 <- kidneyGeneLists[[1]]
#'  gl2 <- kidneyGeneLists[[2]]
#'  orgPkg ='org.Hs.eg.db'
#'  GOIDSinLevel.1 <- GOIDsInLevel (GOLev=GOLev, onto = anOnto,  restricted=FALSE);
#'  length(GOIDSinLevel.1)
#'  GOIDSinLevel.2 <- GOIDsInLevel (GOLev=GOLev, onto = anOnto,  geneList1 = gl1, geneList2 = gl2, orgPackage =orgPkg)
#'  length(GOIDSinLevel.2)
GOIDsInLevel <- function (GOLev, onto, geneList1 = NULL, geneList2 = NULL, orgPackage = NULL,
                          restricted = TRUE){
  levelIDs <- getGOLevel(onto = onto, level = GOLev)
  if (restricted && (!is.null(geneList1)&&(!is.null(geneList2))&&(!is.null(orgPackage)))){
    specGOIDs1 <- GOTermsList(geneList1, onto = onto, orgPkg = orgPackage)
    specGOIDs2 <- GOTermsList(geneList2, onto = onto, orgPkg = orgPackage)
    allGOIDs1 <- unique(unlist(getAncestorsLst(specGOIDs1, onto = onto)))
    allGOIDs2 <- unique(unlist(getAncestorsLst(specGOIDs2, onto = onto)))
    levelIDs <- intersect(levelIDs, union(allGOIDs1, allGOIDs2))
  }
  return(levelIDs)
}

#' crossTabGOIDs
#'
#' This function performs a crosstabulation between two lists of enriched GOTerms
#' The lists are intended to have been obtained from enrichment analyses performed on two gene lists
#'
#' @param GO1 character vector containing a FIRST list of GO identifiers
#' @param GO2 character vector containing a SECOND gene list GO identifiers
#' @param onto string describing the ontology. Belongs to c('BP', 'MF', 'CC', 'ANY')
#' @param GOLev An integer
#' @param listNames character vector with names of the genelists that generated the enriched GOIDs
#' @param geneList1 character vector containing a FIRST gene list of entrez IDs
#' @param geneList2 character vector containing a SECOND gene list of entrez IDs
#' @param orgPackage A string wih the name of the annotation package
#' @param restricted Boolean variable to decide how tabulation is performed.
#' Unrestricted tabulation crosses _all_ GO Terms located at the level indicated by `GOLev` with the two GOIDs lists
#' Restricted tabulation crosses only terms from the selected GO level that are _common to ancestor terms of either list_.
#' That is, if one term in the selected GO level is not an ancestor of at least one of the gene list most specific GO terms
#' it is excluded from the GO Level's terms because it is impossible that it appears as being enriched.
#' @examples
#'  data(kidneyEnrichedGOIDs)
#'  GOIDs1 <-kidneyEnrichedGOIDs[[1]]
#'  GOIDs2 <-kidneyEnrichedGOIDs[[2]]
#'  names4lists<- names(kidneyEnrichedGOIDs)[1:2]
#'  anOnto <- 'BP'
#'  GOLev <- 3
#'  gl1 <- kidneyGeneLists[[1]]
#'  gl2 <- kidneyGeneLists[[2]]
#'  orgPkg ='org.Hs.eg.db'
#'  restrictTab <- TRUE
#'  crossTabbedGOIdsRestricted <- crossTabGOIDs (GO1 = GOIDs1, GO2 = GOIDs2, onto = anOnto,
#'  GOLevel =GOLev, listNames=names4lists, geneList1 = gl1, geneList2 = gl2, orgPackage =orgPkg, restricted = restrictTab)
#'  show(crossTabbedGOIdsRestricted)
#'  restrictTab <- FALSE
#'  crossTabbedGOIdsUnrestricted <- crossTabGOIDs (GO1 = GOIDs1, GO2 = GOIDs2, onto = anOnto,
#'  GOLevel =GOLev, listNames=names4lists, geneList1 = gl1, geneList2 = gl2, orgPackage =orgPkg, restricted = restrictTab)
#'  show(crossTabbedGOIdsUnrestricted)
crossTabGOIDs <- function (GO1, GO2, onto, GOLevel, listNames = NULL,
                           geneList1 = NULL, geneList2 = NULL, orgPackage = NULL,
                           restricted = FALSE)
{
  levelIDs <- GOIDsInLevel (GOLev = GOLevel, onto = onto,  geneList1 = geneList1, geneList2 = geneList2,
                            orgPackage = orgPackage, restricted = restricted)
  levelInGO1 <- levelIDs %in% GO1
  levelInGO2 <- levelIDs %in% GO2
  if (is.null(listNames)) {
    dnnNames <- paste0("Enriched_in_", c("list_1", "list_2"))
  }else{
    dnnNames <- paste0("Enriched_in_", listNames)
  }
  crossedTable <- table(levelInGO1, levelInGO2, dnn=dnnNames)
}

#' crossTab4GeneLists
#'
#' This function builds a cross-tabulation of enriched and non-enriched GO terms
#' from two gene lists
#'
#' @param geneList1 character vector containing a FIRST gene list of entrez IDs
#' @param geneList2 character vector containing a SECOND gene list of entrez IDs
#' @param geneUniverse character vector containing all genes from where geneLists have been extracted
#' @param orgPackage A string wih the name of the annotation package
#' @param onto string describing the ontology. Belongs to c('BP', 'MF', 'CC', 'ANY')
#' @param GOLev An integer
#' @param restricted Boolean variable to decide how tabulation of GOIDs is performed.
#' Unrestricted tabulation crosses _all_ GO Terms located at the level indicated by `GOLev` with the two GOIDs lists
#' Restricted tabulation crosses only terms from the selected GO level that are _common to ancestor terms of either list_.
#' That is, if one term in the selected GO level is not an ancestor of at least one of the gene list most specific GO terms
#' it is excluded from the GO Level's terms because it is impossible that it appears as being enriched.
#' @param pAdjustMeth string describing the adjust method. Belongs to c('BH', 'BY', 'Bonf')
#' @param pvalCutoff A numeric value
#' @param qvalCutoff A numeric value
#' @examples
#'  data(kidneyGeneLists)
#'  data(humanEntrezIDs)
#'  gl1 <-kidneyGeneLists[[1]]
#'  gl2 <-kidneyGeneLists[[2]]
#'  anOnto <- 'BP'
#'  GOLev<- 3
#'  restricted <- FALSE
#'  adjMeth<- 'BH'
#'  pValCut <- 0.05
#'  qValCut <- 0.01
#'  crossTabFromGeneListsUnrestricted <- crossTabGOIDs4GeneLists (genelist1=gl1, genelist2=gl2,
#'  geneUniverse=humanEntrezIDs, orgPackg="org.Hs.eg.db", onto=anOnto, GOLevel=GOLev, restricted=restricted)
#'  restricted <- TRUE
#'  crossTabFromGeneListsRestricted <- crossTabGOIDs4GeneLists (genelist1=gl1, genelist2=gl2,
#'  geneUniverse=humanEntrezIDs, orgPackg="org.Hs.eg.db", onto=anOnto, GOLevel=GOLev, restricted=restricted)
crossTabGOIDs4GeneLists <- function (genelist1, genelist2, geneUniverse, orgPackg,
                                     onto, GOLevel, restricted = FALSE,
                                     pAdjustMeth = "BH", pvalCutoff = 0.01, qvalCutoff = 0.05)
{
  enriched1 <- enrichGO(gene = genelist1, universe = geneUniverse, OrgDb = orgPackg,
                        ont = onto, pvalueCutoff = pvalCutoff)
  GOIDs1 <- as.character(as.data.frame(enriched1)$ID)
  enriched2 <- enrichGO(gene = genelist2, universe = geneUniverse, OrgDb = orgPackg,
                        ont = onto, pvalueCutoff = pvalCutoff)
  GOIDs2 <- as.character(as.data.frame(enriched2)$ID)
  crossTabGOIDs4GeneLists <- crossTabGOIDs (GO1 = GOIDs1, GO2 = GOIDs2,
                                            onto = onto, GOLevel = GOLevel,
                                            geneList1 = genelist1, geneList2 = genelist2,
                                            orgPackage = orgPackg,
                                            restricted = restricted)
}
