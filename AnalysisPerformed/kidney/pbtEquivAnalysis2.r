# ******************************************************************************
#             PBT GENE LISTS JOINT ENRICHMENT EQUIVALENCE ANALYSIS
# ******************************************************************************

load("pbtGeneLists2.rda") # kidney gene lists
library(goSorensen)
library(equivStandardTest)
data(humanEntrezIDs)
source('adjSignifPvals.r')

sapply(pbtGeneLists2, length)

numLists <- length(pbtGeneLists2)
lstNams <- names(pbtGeneLists2)
for (i in 2:length(pbtGeneLists2)) {
  for (j in 1:(i-1)) {
    cat(lstNams[i], "&", lstNams[j], length(intersect(pbtGeneLists2[[i]], pbtGeneLists2[[j]])),
        "common genes of", length(pbtGeneLists2[[i]]), length(pbtGeneLists2[[j]]), "\n")
  }
}

# Some illustrative but preliminary analyzes using package goSorensen:

# Equivalence test at level 2 of the BP ontology between gene lists KT1 and IRITD5
KT1_IRITD5.BP.4 <- equivTestSorensen(pbtGeneLists2[["KT1"]], pbtGeneLists2[["IRITD5"]],
                                     geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
                                     onto = "BP", GOLevel = 4, listNames = c("KT1", "IRITD5"))
names(KT1_IRITD5.BP.4)
KT1_IRITD5.BP.4$statistic
KT1_IRITD5.BP.4$p.value
getPvalue(KT1_IRITD5.BP.4)
KT1_IRITD5.BP.4$conf.int
getUpper(KT1_IRITD5.BP.4)
KT1_IRITD5.BP.4$estimate
KT1_IRITD5.BP.4$null.value # d0
KT1_IRITD5.BP.4$stderr
getSE(KT1_IRITD5.BP.4)
KT1_IRITD5.BP.4$alternative
KT1_IRITD5.BP.4$method
KT1_IRITD5.BP.4$enrichTab
getTable(KT1_IRITD5.BP.4)

# Number of annotated GO items in the union of both gene sets, in ontology BP at level 4:
s <- sum(KT1_IRITD5.BP.4$enrichTab)
names(s) <- "annotated GO items"
s

# All pairwise equivalence tests at level 5 of the BP ontology
BP.5 <- equivTestSorensen(pbtGeneLists2,
                          geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
                          onto = "BP", GOLevel = 5)
BP.5


# All pairwise equivalence tests for all levels from 3 to 10 and for all GO ontologies
# (Extremely time consuming. To save time jump, uncomment and run 'load("pbtEquivSorensen2.rda")'

pbtAllOntosAndLevels2 <- allEquivTestSorensen(pbtGeneLists2,
                                              geneUniverse = humanEntrezIDs,
                                              orgPackg = "org.Hs.eg.db")
# save(pbtAllOntosAndLevels2, file = "pbtEquivSorensen2.rda")
# load("pbtEquivSorensen2.rda")

# From 91 possible p-values (91 = 14 * (14 - 1) / 2) and after the Holm's adjustment for testing multiplicity,
# identify those who are <= 0.05?, excluding NA values:
# ---------------------------------------------------------------------------------------------------------
# For ontology BP:
# Under d0 = 0.4444
sapply(getPvalue(pbtAllOntosAndLevels2, onto = "BP")[["BP"]], function(thisLevPvals){
  thisLevPvals <- p.adjust(thisLevPvals[!is.na(thisLevPvals)], method = "holm")
  thisLevPvals[thisLevPvals <= 0.05]
})
# Under a more restrictive d0 = 0.2857
sapply(getPvalue(upgrade(pbtAllOntosAndLevels2, d0 = 1/(1 + 2*1.25)), onto = "BP")[["BP"]], function(thisLevPvals){
  thisLevPvals <- p.adjust(thisLevPvals[!is.na(thisLevPvals)], method = "holm")
  thisLevPvals[thisLevPvals <= 0.05]
})

# CAUTION! Some of these "significant" results have a very low reliability: the enrichment frequencies
# are extremely low:
getTable(pbtAllOntosAndLevels2, onto = "BP")

# ---------------------------------------------------------------------------------------------------------
# For ontology CC:
# Under d0 = 0.4444
sapply(getPvalue(pbtAllOntosAndLevels2, onto = "CC")[["CC"]], function(thisLevPvals){
  thisLevPvals <- p.adjust(thisLevPvals[!is.na(thisLevPvals)], method = "holm")
  thisLevPvals[thisLevPvals <= 0.05]
})
# Under a more restrictive d0 = 0.2857
sapply(getPvalue(upgrade(pbtAllOntosAndLevels2, d0 = 1/(1 + 2*1.25)), onto = "CC")[["CC"]], function(thisLevPvals){
  thisLevPvals <- p.adjust(thisLevPvals[!is.na(thisLevPvals)], method = "holm")
  thisLevPvals[thisLevPvals <= 0.05]
})

# CAUTION! Some of these "significant" results have a very low reliability: the enrichment frequencies
# are extremely low:
getTable(pbtAllOntosAndLevels2, onto = "CC")


# ---------------------------------------------------------------------------------------------------------
# For ontology MF:
# Under d0 = 0.4444
sapply(getPvalue(pbtAllOntosAndLevels2, onto = "MF")[["MF"]], function(thisLevPvals){
  thisLevPvals <- p.adjust(thisLevPvals[!is.na(thisLevPvals)], method = "holm")
  thisLevPvals[thisLevPvals <= 0.05]
})
# Under a more restrictive d0 = 0.2857
sapply(getPvalue(upgrade(pbtAllOntosAndLevels2, d0 = 1/(1 + 2*1.25)), onto = "MF")[["MF"]], function(thisLevPvals){
  thisLevPvals <- p.adjust(thisLevPvals[!is.na(thisLevPvals)], method = "holm")
  thisLevPvals[thisLevPvals <= 0.05]
})

# CAUTION! Some of these "significant" results have a very low reliability: the enrichment frequencies
# are extremely low:
getTable(pbtAllOntosAndLevels2, onto = "MF")

# CAUTION! Some of these "significant" results may have a very low reliability if the joint enrichment
# frequencies are extremely low.
# This function returns the adjusted significant p-values jointly with the enrichment contingency tables,
# in order to put these values in an adequate context (e.g., those who are not very credible):
signifPvals_d0_0.4444 <- adjSignifPvals(pbtAllOntosAndLevels2)


# In BP
signifPvals_d0_0.4444$BP
# In CC,
signifPvals_d0_0.4444$CC
# In MF,
signifPvals_d0_0.4444$MF
sink(file = NULL)
# For a more restrictive d0 = 0.2857:
signifPvals_d0_0.2857 <- adjSignifPvals(upgrade(cancerEquivSorensen, d0 = 1/(1 + 2*1.25)))

# In BP ontology,
signifPvals_d0_0.2857$BP
# In CC,
signifPvals_d0_0.2857$CC
# In MF:
signifPvals_d0_0.2857$MF
