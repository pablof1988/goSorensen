## ----style, echo = FALSE, results = 'asis'------------------------------------
BiocStyle::markdown()

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(dpi=100,fig.width=7)

## ----env, message = FALSE, warning = FALSE, echo = TRUE-----------------------
library(goSorensen)

## -----------------------------------------------------------------------------
data("pbtGeneLists")

## ----comment=NA---------------------------------------------------------------
sapply(pbtGeneLists, length)

## ----warning=FALSE, message=FALSE, comment=NA---------------------------------
# Previously load the genomic annotation package for the studied specie:
library(org.Hs.eg.db)
humanEntrezIDs <- keys(org.Hs.eg.db, keytype = "ENTREZID")

# Computing the irrelevance-threshold matrix of dissimilarities
dismatBP5 <- sorenThreshold(pbtGeneLists, onto = "BP", GOLevel = 5,
                            geneUniverse = humanEntrezIDs, 
                            orgPackg = "org.Hs.eg.db")
round(dismatBP5, 2)

## ----warning=FALSE, message=FALSE, comment=NA, eval=TRUE----------------------
# Previously compute the 2x2 contingency tables:
ctableBP5 <- buildEnrichTable(pbtGeneLists, onto = "BP", GOLevel = 5,
                            geneUniverse = humanEntrezIDs, 
                            orgPackg = "org.Hs.eg.db")

# Computing the irrelevance-threshold matrix of dissimilarities
dismatBP5 <- sorenThreshold(ctableBP5)

## ----warning=FALSE, message=FALSE, comment=NA, eval=FALSE---------------------
# allDismat <- allSorenThreshold (pbtGeneLists, geneUniverse = humanEntrezIDs,
#                             orgPackg = "org.Hs.eg.db", ontos = c("BP", "CC"),
#                             GOLevels = 4:7)

## ----warning=FALSE, message=FALSE, comment=NA, eval=FALSE---------------------
# allTabs <- allBuildEnrichTable(pbtGeneLists, geneUniverse = humanEntrezIDs,
#                             orgPackg = "org.Hs.eg.db", ontos = c("BP", "CC"),
#                             GOLevels = 4:7)
# allDismat <- allSorenThreshold(allTabs)

## ----warning=FALSE, message=FALSE, comment=NA---------------------------------
clust.threshold <- hclustThreshold(dismatBP5)
plot(clust.threshold)

## ----warning=FALSE, message=FALSE, comment=NA---------------------------------
mds <- as.data.frame(cmdscale(dismatBP5, k = 2))
mds

## ----warning=FALSE, message=FALSE, comment=NA, fig.align='center'-------------
library(ggplot2)
library(ggrepel)
graph <- ggplot() +
  geom_point(aes(mds[,1], mds[,2]), color = "blue", size = 3) +
  geom_text_repel(aes(mds[,1], mds[,2], label = attr(dismatBP5, "Labels")), 
                  color = "black", size = 3) +
    xlab("Dim 1") +
    ylab("Dim 2") +
    theme_light()
graph

## ----warning=FALSE, message=FALSE, comment=NA---------------------------------
# Split the axis 20% to the left, 60% to the middle and 20% to the right:
prop <- c(0.2, 0.6, 0.2) 
# Sort according  dimension 1:
sorted <- mds[order(mds[, 1]), ] 
# Determine the range for dimension 1.
range <- sorted[, 1][c(1, nrow(mds))] 
# Find the cutpoints to split the axis:
cutpoints <- cumsum(prop) * diff(range) + range[1]
cutpoints <- cutpoints[-length(cutpoints)]

# Identify lists to the left:  
lleft <- rownames(sorted[sorted[, 1] < cutpoints[1], ])
# Identify lists to the right
lright <- rownames(sorted[sorted[, 1] > cutpoints[2], ]) 

lleft
lright

## ----warning=FALSE, message=FALSE, comment=NA, fig.align='center'-------------
graph +
  geom_vline(xintercept = cutpoints, color = "red", linetype = "dashed", 
             linewidth = 0.75)

## -----------------------------------------------------------------------------
tableleft <- attr(ctableBP5, "enriched")[, lleft]

tableright <- attr(ctableBP5, "enriched")[, lright]

## ----warning=FALSE, message=FALSE, comment=NA---------------------------------
tableleft[1:15, ]
tableright[1:15, ]

## ----warning=FALSE, message=FALSE, comment=NA---------------------------------
lmnsd <- apply(tableleft, 1, 
                 function(x){c("meanLeft" = mean(x), "sdLeft" = sd(x))})
rmnsd <- apply(tableright, 1, 
                 function(x){c("meanRight" = mean(x), "sdRight" = sd(x))})

## ----warning=FALSE, message=FALSE, comment=NA---------------------------------
nl <- ncol(tableleft) 
nr <- ncol(tableright) 
stat <- abs(lmnsd[1, ] - rmnsd[1, ]) / 
  sqrt((((lmnsd[2, ] / nl) + (rmnsd[2, ] / nr))) + 0.00000001)

## ----warning=FALSE, message=FALSE, comment=NA---------------------------------
result <- stat[stat == max(stat)]
result

## ----warning=FALSE, message=FALSE, comment=NA---------------------------------
library(GO.db)
library(knitr)
# Previous function to get the description of the identified GO terms:
get_go_description <- function(go_id) {
  go_term <- Term(GOTERM[[go_id]])
  return(go_term)
}

# GO terms description:
kable(data.frame(Description = sapply(names(result), get_go_description, 
                                             USE.NAMES = TRUE)))

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

