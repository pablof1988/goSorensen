library(goSorensen)
library(equivStandardTest)
source("adjSignifPvals.R")
data(humanEntrezIDs)


# package goSorensen authomatically charges object "allOncoGeneLists"
allOncoGeneLists
sapply(allOncoGeneLists, length)

# Formerly, the original set of gene lists was reduced to those with almost 100 genes:
# allOncoGeneLists <- allOncoGeneLists[sapply(allOncoGeneLists, length) >= 100]
# sapply(allOncoGeneLists, length)

numLists <- length(allOncoGeneLists)
lstNams <- names(allOncoGeneLists)
for (i in 2:length(allOncoGeneLists)) {
  for (j in 1:(i-1)) {
    cat(lstNams[i], "&", lstNams[j], length(intersect(allOncoGeneLists[[i]], allOncoGeneLists[[j]])),
        "common genes of", length(allOncoGeneLists[[i]]), length(allOncoGeneLists[[j]]), "\n")
  }
}


# The generic function "equivTestSorensen" in package "goSorensen" implements the equivalence test based
# on the Sorensen-Dice distance.
# There are many methods of this function, for different classes of objects passed as arguments.
# Providing two gene lists (essentially, to "character" vectors of gene identifiers) as arguments, it returns
# an object of class "equivSDhtest" inheriting from "htest".
# Providing an object of class "list" of length k with each of its elements representing a gene list (i.e., a "list"
# of k "character" vectors), all possible k*(k - 1)/2 pairwise tests are performed and an object of class "equivSDhtestList"
# (in fact, a "list") of "htest" objects is generated.
# There are also methods for data already summarized in form of 2x2 contingency tables of joint enrichment.

# Examples:
# Equivalence test between gene lists 'waldman' and 'atlas', in dataset 'cancerGeneLists',
# at level 4 of the BP ontology:
# (Slightly time consuming, you can jump this statement because package 'goSorensen' authomatically
# charges 'waldman_atlas.BP.4' dataset)
waldman_atlas.BP.4 <- equivTestSorensen(allOncoGeneLists[["waldman"]], allOncoGeneLists[["atlas"]],
                                        geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
                                        onto = "BP", GOLevel = 4, listNames = c("waldman", "atlas"))
waldman_atlas.BP.4
class(waldman_atlas.BP.4)
names(waldman_atlas.BP.4)

waldman_atlas.BP.4$statistic
waldman_atlas.BP.4$p.value
# or:
getPvalue(waldman_atlas.BP.4)

waldman_atlas.BP.4$estimate
# or
getDissimilarity(waldman_atlas.BP.4)

waldman_atlas.BP.4$stderr
# or:
getSE(waldman_atlas.BP.4)

waldman_atlas.BP.4$conf.int
# or:
getUpper(waldman_atlas.BP.4)

getDissimilarity(waldman_atlas.BP.4) + qnorm(0.95) * getSE(waldman_atlas.BP.4)

waldman_atlas.BP.4$null.value # d0

waldman_atlas.BP.4$alternative
waldman_atlas.BP.4$method

waldman_atlas.BP.4$enrichTab
# or:
getTable(waldman_atlas.BP.4)

upgrade(waldman_atlas.BP.4, d0 = 1/(1 + 10/9)) # d0 = 0.4737
waldman_atlas.BP.4_strict <- upgrade(waldman_atlas.BP.4, d0 = 1/(1 + 2*1.25)) # d0 = 0.2857
waldman_atlas.BP.4_strict
class(waldman_atlas.BP.4_strict)
upgrade(waldman_atlas.BP.4, d0 = 1/(1 + 2*1.25), conf.level = 0.99)

# Number of annotated GO items in the union of both gene sets, in ontology BP at level 4:
s <- sum(getTable(waldman_atlas.BP.4))
names(s) <- "annotated GO items"
s

# All pairwise equivalence tests at level 4 of the BP ontology (quite time consuming, you can jump
# this sentence and use the dataset 'BP.4' which is directly charged with package 'goSorensen')
BP.4 <- equivTestSorensen(allOncoGeneLists,
                          geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
                          onto = "BP", GOLevel = 4)

class(BP.4)
BP.4
getPvalue(BP.4)
getPvalue(BP.4, simplify = FALSE)

getDissimilarity(BP.4)
getDissimilarity(BP.4, simplify = FALSE)

getUpper(BP.4)
getUpper(BP.4, simplify = FALSE)

getDissimilarity(BP.4, simplify = FALSE) + qnorm(0.95) * getSE(BP.4, simplify = FALSE)

getSE(BP.4)
getSE(BP.4, simplify = FALSE)
getTable(BP.4)

BP.4_strict <- upgrade(BP.4, d0 = 1/(1 + 2*1.25)) # d0 = 0.2857
BP.4_strict
class(BP.4_strict)

upgrade(BP.4, d0 = 1/(1 + 2*1.25),
        conf.level = 0.90, check.table = FALSE)

# (Very time consuming. Alternatively, you may use the dataset 'cancerEquivSorensen' directly,
# it is automatically charged with the package 'goSorensen'),
# By default, the tests are iterated over all GO ontologies and for levels 3 to 10: (these
# results are available at  cancerEquivSorensen.rda file)
cancerEquivSorensen <- allEquivTestSorensen(allOncoGeneLists,
                                            geneUniverse = humanEntrezIDs,
                                            orgPackg = "org.Hs.eg.db")
#save(cancerEquivSorensen, file = "cancerEquivSorensen.rda")
#load('cancerEquivSorensen.rda')
# # (Also very time consuming.) It is not required to iterate over all ontologies,
# # "allEquivTestSorensen" iterates these procedures over the specified GO ontologies and levels,
# # e.g., to iterate only over the GO ontologies MF and BP and levels 5 and 6:
# allEquivTestSorensen(allOncoGeneLists, ontos = c("MF", "BP"), GOLevels = 5:6,
#                      geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
# # With a more strict equivalence limit d0:
# allEquivTestSorensen(allOncoGeneLists, ontos = c("MF", "BP"), GOLevels = 5:6, d0 = 1 / (1 + 2 * (10/9)), #d0 =0.3103
#                      geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")

cancerEquivSorensen

# From 21 possible p-values (21 = 7 * (7 - 1) / 2) and after the Holm's adjustment for testing multiplicity,
# identify those who are <= 0.05? (Excluding NA values)
#   BUT JUMP TO LINE 248 FOR A MORE INFORMATIVE OUTPUT, INCLUDING THE 2x2 CONTINGENCY TABLES OF ENRICHMENT
# ---------------------------------------------------------------------------------------------------------
# For ontology BP:
# Under d0 = 0.4444
sapply(getPvalue(cancerEquivSorensen, onto = "BP")[["BP"]], function(thisLevPvals){
  thisLevPvals <- p.adjust(thisLevPvals[!is.na(thisLevPvals)], method = "holm")
  thisLevPvals[thisLevPvals <= 0.05]
})

# Under a more restrictive d0 = 0.2857
sapply(getPvalue(upgrade(cancerEquivSorensen, d0 = 1/(1 + 2*1.25)), onto = "BP")[["BP"]], function(thisLevPvals){
  thisLevPvals <- p.adjust(thisLevPvals[!is.na(thisLevPvals)], method = "holm")
  thisLevPvals[thisLevPvals <= 0.05]
})
# ************** Highly stable results along all GO levels for BP ontology **************************
# ************** a set of equivalent lists seems to be identified          **************************

# ---------------------------------------------------------------------------------------------------------
# For ontology CC:
# Under d0 = 0.4444
sapply(getPvalue(cancerEquivSorensen, onto = "CC")[["CC"]], function(thisLevPvals){
  thisLevPvals <- p.adjust(thisLevPvals[!is.na(thisLevPvals)], method = "holm")
  thisLevPvals[thisLevPvals <= 0.05]
})

# Under a more restrictive d0 = 0.2857
sapply(getPvalue(upgrade(cancerEquivSorensen, d0 = 1/(1 + 2*1.25)), onto = "CC")[["CC"]], function(thisLevPvals){
  thisLevPvals <- p.adjust(thisLevPvals[!is.na(thisLevPvals)], method = "holm")
  thisLevPvals[thisLevPvals <= 0.05]
})
# **** Stable but possibly less interesting (less similarity between lists) results for CC ontology ****

# ---------------------------------------------------------------------------------------------------------
# For ontology MF:
# Under d0 = 0.4444
sapply(getPvalue(cancerEquivSorensen, onto = "MF")[["MF"]], function(thisLevPvals){
  thisLevPvals <- p.adjust(thisLevPvals[!is.na(thisLevPvals)], method = "holm")
  thisLevPvals[thisLevPvals <= 0.05]
})

# Under a more restrictive d0 = 0.2857
sapply(getPvalue(upgrade(cancerEquivSorensen, d0 = 1/(1 + 2*1.25)), onto = "MF")[["MF"]], function(thisLevPvals){
  thisLevPvals <- p.adjust(thisLevPvals[!is.na(thisLevPvals)], method = "holm")
  thisLevPvals[thisLevPvals <= 0.05]
})

# **** Stable but possibly less interesting (less similarity between lists) results for MF ontology ****
# ---------------------------------------------------------------------------------------------------------

class(cancerEquivSorensen$BP$`level 4`)
class(cancerEquivSorensen$BP$`level 4`$humanlymph)
class(cancerEquivSorensen$BP$`level 4`$humanlymph$atlas)

# 2x2 contingecy tables of joint enrichment:
getTable(cancerEquivSorensen)
getTable(cancerEquivSorensen, GOLevel = "level 6")
getTable(cancerEquivSorensen, GOLevel = "level 6", listNames = c("waldman", "sanger"))
getTable(cancerEquivSorensen, GOLevel = "level 6", onto = "BP")
getTable(cancerEquivSorensen, GOLevel = "level 6", onto = "BP", listNames = c("waldman", "sanger"))
getTable(cancerEquivSorensen$BP$`level 4`)

# p-values:
getPvalue(cancerEquivSorensen)
getPvalue(cancerEquivSorensen, simplify = FALSE)
getPvalue(cancerEquivSorensen, GOLevel = "level 6")
getPvalue(cancerEquivSorensen, GOLevel = "level 6", simplify = FALSE)
getPvalue(cancerEquivSorensen, GOLevel = "level 6", listNames = c("waldman", "sanger"))
getPvalue(cancerEquivSorensen, GOLevel = "level 6", onto = "BP")
getPvalue(cancerEquivSorensen, GOLevel = "level 6", onto = "BP", simplify = FALSE)
getPvalue(cancerEquivSorensen, GOLevel = "level 6", onto = "BP", listNames = c("waldman", "sanger"))
getPvalue(cancerEquivSorensen$BP$`level 4`)

# Sorensen-Dice dissimilarity:
getDissimilarity(cancerEquivSorensen)
getDissimilarity(cancerEquivSorensen, simplify = FALSE)
getDissimilarity(cancerEquivSorensen, GOLevel = "level 6")
getDissimilarity(cancerEquivSorensen, GOLevel = "level 6", simplify = FALSE)
getDissimilarity(cancerEquivSorensen, GOLevel = "level 6", listNames = c("waldman", "sanger"))
getDissimilarity(cancerEquivSorensen, GOLevel = "level 6", onto = "BP")
getDissimilarity(cancerEquivSorensen, GOLevel = "level 6", onto = "BP", simplify = FALSE)
getDissimilarity(cancerEquivSorensen, GOLevel = "level 6", onto = "BP", listNames = c("waldman", "sanger"))
getDissimilarity(cancerEquivSorensen$BP$`level 4`)

# Upper confidence limits:
getUpper(cancerEquivSorensen)
getUpper(cancerEquivSorensen, simplify = FALSE)
getUpper(cancerEquivSorensen, GOLevel = "level 6")
getUpper(cancerEquivSorensen, GOLevel = "level 6", simplify = FALSE)
getUpper(cancerEquivSorensen, GOLevel = "level 6", listNames = c("waldman", "sanger"))
getUpper(cancerEquivSorensen, GOLevel = "level 6", onto = "BP")
getUpper(cancerEquivSorensen, GOLevel = "level 6", onto = "BP", simplify = FALSE)
getUpper(cancerEquivSorensen, GOLevel = "level 6", onto = "BP", listNames = c("waldman", "sanger"))
getUpper(cancerEquivSorensen$BP$`level 4`)

# Standard error of the Sorensen-Dice dissimilarity estimate:
getSE(cancerEquivSorensen)
getSE(cancerEquivSorensen, simplify = FALSE)
getSE(cancerEquivSorensen, GOLevel = "level 6")
getSE(cancerEquivSorensen, GOLevel = "level 6", simplify = FALSE)
getSE(cancerEquivSorensen, GOLevel = "level 6", listNames = c("waldman", "sanger"))
getSE(cancerEquivSorensen, GOLevel = "level 6", onto = "BP")
getSE(cancerEquivSorensen, GOLevel = "level 6", onto = "BP", simplify = FALSE)
getSE(cancerEquivSorensen, GOLevel = "level 6", onto = "BP", listNames = c("waldman", "sanger"))
getSE(cancerEquivSorensen$BP$`level 4`)

cancerEquivSorensen2 <- upgrade(cancerEquivSorensen, d0 = 1/(1 + 2*1.25)) # d0 = 0.2857
cancerEquivSorensen2
class(cancerEquivSorensen2)
class(cancerEquivSorensen2$BP$`level 4`)
class(cancerEquivSorensen$BP$`level 4`$humanlymph)
class(cancerEquivSorensen$BP$`level 4`$humanlymph$atlas)

# CAUTION! Some of these "significant" results may have a very low reliability if the joint enrichment
# frequencies are extremely low.
# This function returns the adjusted significant p-values jointly with the enrichment contingency tables,
# in order to put these values in an adequate context (e.g., those who are not very credible):
signifPvals_d0_0.4444 <- adjSignifPvals(cancerEquivSorensen)

# In BP ontology, all detected equivalencies seem very reliable:
signifPvals_d0_0.4444$BP
# In CC, this is not so clear. The few significant equivalencies are associated to low enrichment
# frequencies:
signifPvals_d0_0.4444$CC
# In MF, the constantly present equivalence between Vogelstein and sanger seems especially reliable
# only at levels 4 and 5:
signifPvals_d0_0.4444$MF

# For a more restrictive d0 = 0.2857:
signifPvals_d0_0.2857 <- adjSignifPvals(upgrade(cancerEquivSorensen, d0 = 1/(1 + 2*1.25)))

# In BP ontology, a consistent and reliable set of equivalences is still present, at all GO levels:
signifPvals_d0_0.2857$BP
# In CC, even less significant equivalencies and not very reliable, only for the first GO levels:
signifPvals_d0_0.2857$CC
# In MF, again the same results: a persistent equivalence between Vogelstein and sanger, at GO
# levels 4 and 5:
signifPvals_d0_0.2857$MF
