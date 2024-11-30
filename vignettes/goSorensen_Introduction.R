## ----style, echo = FALSE, results = 'asis'------------------------------------
BiocStyle::markdown()

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(dpi=25,fig.width=7)

## ----env, message = FALSE, warning = FALSE, echo = TRUE-----------------------
library(goSorensen)

## ---- eval=TRUE---------------------------------------------------------------
if (!requireNamespace("goSorensen", quietly = TRUE)) {
    BiocManager::install("goSorensen")
}
library(goSorensen)

## -----------------------------------------------------------------------------
data("allOncoGeneLists")

## -----------------------------------------------------------------------------
length(allOncoGeneLists)
sapply(allOncoGeneLists, length)

# First 20 gene identifiers of gene lists Vogelstein and sanger:
allOncoGeneLists[["Vogelstein"]][1:20]
allOncoGeneLists[["sanger"]][1:20]

## ---- message = FALSE, warning = FALSE, eval = TRUE---------------------------
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    BiocManager::install("org.Hs.eg.db")
}
library(org.Hs.eg.db)

## ---- message = FALSE, warning = FALSE, echo = TRUE---------------------------
humanEntrezIDs <- keys(org.Hs.eg.db, keytype = "ENTREZID")

## ---- warning=FALSE, message=FALSE--------------------------------------------
# Build the enrichment contingency table between gene lists Vogelstein and 
# sanger for the MF ontology at GO level 5:
enrichTab <- buildEnrichTable(allOncoGeneLists[["Vogelstein"]],
                              allOncoGeneLists[["sanger"]],
                              geneUniverse = humanEntrezIDs, 
                              orgPackg = "org.Hs.eg.db",
                              onto = "MF", GOLevel = 5, 
                              listNames = c("Vogelstein", "sanger"))
enrichTab

# Equivalence test for an equivalence (or negligibility) limit 0.2857
testResult <- equivTestSorensen(enrichTab, d0 = 0.2857)
testResult

## -----------------------------------------------------------------------------
equivTestSorensen(allOncoGeneLists[["Vogelstein"]],
                  allOncoGeneLists[["sanger"]], d0 = 0.2857,
                  geneUniverse = humanEntrezIDs, 
                  orgPackg = "org.Hs.eg.db",
                  onto = "MF", GOLevel = 5, 
                  listNames = c("Vogelstein", "sanger"))

## -----------------------------------------------------------------------------
boot.testResult <- equivTestSorensen(enrichTab, d0 = 0.2857, boot = TRUE)
boot.testResult

## -----------------------------------------------------------------------------
# The Sorensen dissimilarity from the contingency table:
dSorensen(enrichTab)
# The Sorensen dissimilarity from the gene lists:
dSorensen(allOncoGeneLists[["Vogelstein"]], 
          allOncoGeneLists[["sanger"]],
          geneUniverse = humanEntrezIDs, 
          orgPackg = "org.Hs.eg.db",
          onto = "MF", GOLevel = 5, 
          listNames = c("Vogelstein", "sanger"))

# The standard error from the contingency table::
seSorensen(enrichTab)
# or from the gene lists:
seSorensen(allOncoGeneLists[["Vogelstein"]], 
           allOncoGeneLists[["sanger"]],
           geneUniverse = humanEntrezIDs, 
           orgPackg = "org.Hs.eg.db",
           onto = "MF", GOLevel = 5, 
           listNames = c("Vogelstein", "sanger"))

# Upper limit of the confidence interval from the contingency table:
duppSorensen(enrichTab)
duppSorensen(enrichTab, conf.level = 0.90)
duppSorensen(enrichTab, conf.level = 0.90, boot = TRUE)

# Upper limit of the confidence interval from the gene lists:
duppSorensen(allOncoGeneLists[["Vogelstein"]], 
             allOncoGeneLists[["sanger"]],    
             geneUniverse = humanEntrezIDs, 
             orgPackg = "org.Hs.eg.db", 
             onto = "MF", GOLevel = 5, 
             listNames = c("Vogelstein", "sanger"))

## -----------------------------------------------------------------------------
getDissimilarity(testResult)
getSE(testResult)
getPvalue(testResult)
getTable(testResult)
getUpper(testResult)

# In the bootstrap approach, only these differ:
getPvalue(boot.testResult)
getUpper(boot.testResult)
# (Only available for bootstrap tests) efective number of bootstrap resamples:
getNboot(boot.testResult)

## -----------------------------------------------------------------------------
upgrade(testResult, d0 = 0.4444, conf.level = 0.99, boot = TRUE)

## -----------------------------------------------------------------------------
totalDiss <- dSorensen(allOncoGeneLists, onto = "MF", GOLevel = 5, 
          geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
round(totalDiss, 2)

## -----------------------------------------------------------------------------
allTests <- equivTestSorensen(allOncoGeneLists, d0 = 0.2857, 
                              onto = "MF", GOLevel = 5, 
                              geneUniverse = humanEntrezIDs, 
                              orgPackg = "org.Hs.eg.db")
getPvalue(allTests)
getDissimilarity(allTests, simplify = FALSE)

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

