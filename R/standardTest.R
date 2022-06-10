#' @importFrom clusterProfiler enrichGO
#' @importFrom goProfiles getGOLevel GOTermsList getAncestorsLst
#' @importFrom stats fisher.test
#' @import GO.db
#' @import org.Hs.eg.db

extractGOIDs <- function (enriched) {
  as.character(as.data.frame(enriched)$ID)
}

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
#
# stdTest4GOIDs <- function (GO1, GO2, onto, GOLevel, listNames=NULL)
# {
#   crossTabbedGOIds <- crossTabGOIDs (GO1 = GOIDs1, GO2 = GOIDs2, onto = onto,
#                                      GOLevel =GOLev, listNames=names4lists, restricted=FALSE)
#   fisher.test (crossTabbedGOIds, alt="g")
# }

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

stdTest4GeneLists <- function (genelist1, genelist2, geneUniverse, orgPackg, onto, GOLevel, 
                               restricted = FALSE,
                               pAdjustMeth = "BH", pvalCutoff = 0.01, qvalCutoff = 0.05)
{
  enriched1 <- enrichOnto(geneList = genelist1, geneUniverse = geneUniverse,
                          orgPackage = orgPackg, onto = onto,
                          pvalCutoff = pvalCutoff, qvalCutoff = qvalCutoff, pAdjustMeth = pAdjustMeth)
  GOIDs1 <- as.character(as.data.frame(enriched1)$ID)
  enriched2 <- enrichOnto(geneList = genelist2, geneUniverse = geneUniverse, orgPackage = orgPackg, 
                          onto = onto,
                          pvalCutoff = pvalCutoff, qvalCutoff = qvalCutoff, pAdjustMeth = pAdjustMeth)
  GOIDs2 <- as.character(as.data.frame(enriched2)$ID)
  crossTabbedGOIDs4GeneLists <- crossTabGOIDs (GO1 = GOIDs1, GO2 = GOIDs2, 
                                               onto = onto, GOLevel = GOLevel,
                                               geneList1 = genelist1, geneList2 = genelist2, 
                                               orgPackage = orgPackg,
                                               restricted = restricted)
  fisher.test (crossTabbedGOIDs4GeneLists, alternative = "greater")
}

