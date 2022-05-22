library(goSorensen)

# A contingency table of GO terms mutual enrichment
# between gene lists "atlas" and "sanger":
tab_atlas.sanger_BP3
?tab_atlas.sanger_BP3
class(tab_atlas.sanger_BP3)

# Sorensen-Dice dissimilarity on this contingency table:
?dSorensen
dSorensen(tab_atlas.sanger_BP3)

# Standard error of this Sorensen-Dice dissimilarity estimate:
?seSorensen
seSorensen(tab_atlas.sanger_BP3)

# Upper 95% confidence limit for the Sorensen-Dice dissimilarity:
?duppSorensen
duppSorensen(tab_atlas.sanger_BP3)
# This confidence limit is based on an assimptotic normal N(0,1)
# approximation to the distribution of (dSampl - d) / se, where
# dSampl stands for the sample dissimilarity, d for the true dissimilarity
# and se for the sample dissimilarity standard error estimate.

# Upper confidence limit but using a Student's t instead of a N(0,1)
# (just as an example, not recommended -no theoretical justification)
df <- sum(tab_atlas.sanger_BP3[1:3]) - 2
duppSorensen(tab_atlas.sanger_BP3, z.conf.level = qt(1 - 0.95, df))

# Upper confidence limit but using a bootstrap approximation
# to the sampling distribution, instead of a N(0,1)
set.seed(123)
duppSorensen(tab_atlas.sanger_BP3, boot = TRUE)


# Some computations on diverse data structures:
badConti <- as.table(matrix(c(501, 27, 36, 12, 43, 15, 0, 0, 0),
                            nrow = 3, ncol = 3,
                            dimnames = list(c("a1","a2","a3"),
                                            c("b1", "b2","b3"))))
nice2x2Table(badConti)
incompleteConti <- badConti[1,1:min(2,ncol(badConti)), drop = FALSE]
incompleteConti
nice2x2Table(incompleteConti)

contiAsVector <- c(32, 21, 81, 1439)
nice2x2Table(contiAsVector)
contiAsVector.mat <- matrix(contiAsVector, nrow = 2)
contiAsVector.mat
contiAsVectorLen3 <- c(32, 21, 81)
nice2x2Table(contiAsVectorLen3)

try(dSorensen(badConti), TRUE)
# Apparently it works fine, but returns a wrong value!
dSorensen(badConti, check.table = FALSE)

try(dSorensen(incompleteConti), TRUE)
dSorensen(contiAsVector)
dSorensen(contiAsVector.mat)
dSorensen(contiAsVectorLen3)
dSorensen(contiAsVectorLen3, check.table = FALSE)

try(seSorensen(badConti), TRUE)
try(seSorensen(incompleteConti), TRUE)
seSorensen(contiAsVector)
seSorensen(contiAsVector.mat)
seSorensen(contiAsVectorLen3)
seSorensen(contiAsVectorLen3, check.table = FALSE)

try(duppSorensen(badConti), TRUE)
try(duppSorensen(incompleteConti), TRUE)
duppSorensen(contiAsVector)
duppSorensen(contiAsVector.mat)
set.seed(123)
duppSorensen(contiAsVector, boot = TRUE)
set.seed(123)
duppSorensen(contiAsVector.mat, boot = TRUE)
duppSorensen(contiAsVectorLen3)
# Bootstrapping requires full contingency tables (4 values)
set.seed(123)
try(duppSorensen(contiAsVectorLen3, boot = TRUE), TRUE)

# Equivalence test, H0: d >= d0 vs  H1: d < d0 (d0 = 0.4444)
?equivTestSorensen
equiv.atlas.sanger <- equivTestSorensen(tab_atlas.sanger_BP3)
equiv.atlas.sanger
getTable(equiv.atlas.sanger)
getPvalue(equiv.atlas.sanger)

try(equivTestSorensen(badConti), TRUE)
try(equivTestSorensen(incompleteConti), TRUE)
equivTestSorensen(contiAsVector)
equivTestSorensen(contiAsVector.mat)
set.seed(123)
equivTestSorensen(contiAsVector.mat, boot = TRUE)
equivTestSorensen(contiAsVectorLen3)

try(equivTestSorensen(contiAsVectorLen3, boot = TRUE), TRUE)

# Sorensen-Dice computations from scratch, directly from gene lists
?allOncoGeneLists
# First, the mutual GO node enrichment tables are built, then computations
# proceed from these contingency tables.
# Building the contingency tables is a slow process (many enrichment tests)
normTest <- equivTestSorensen(allOncoGeneLists[["atlas"]], allOncoGeneLists[["sanger"]],
                              listNames = c("atlas", "sanger"),
                              onto = "BP", GOLevel = 5,
                              geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
normTest

# To perform a bootstrap test from scratch would be even slower:
# set.seed(123)
# bootTest <- equivTestSorensen(allOncoGeneLists[["atlas"]], allOncoGeneLists[["sanger"]],
#                               listNames = c("atlas", "sanger"),
#                               boot = TRUE,
#                               onto = "BP", GOLevel = 5,
#                               geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
# bootTest

# It is much faster to upgrade 'normTest' to be a bootstrap test:
set.seed(123)
bootTest <- upgrade(normTest, boot = TRUE)
bootTest
# To know the number of valid bootstrap replicates:
getNboot(bootTest)

# There are similar methods for dSorensen, seSorensen, duppSorensen, etc. to
# compute directly from a pair of gene lists.
# They are quite slow for the same reason as before (many enrichment tests).
# dSorensen(allOncoGeneLists[["atlas"]], allOncoGeneLists[["sanger"]],
#           listNames = c("atlas", "sanger"),
#           onto = "BP", GOLevel = 5,
#           geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
# seSorensen(allOncoGeneLists[["atlas"]], allOncoGeneLists[["sanger"]],
#            listNames = c("atlas", "sanger"),
#            onto = "BP", GOLevel = 5,
#            geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
#
# duppSorensen(allOncoGeneLists[["atlas"]], allOncoGeneLists[["sanger"]],
#              listNames = c("atlas", "sanger"),
#              onto = "BP", GOLevel = 5,
#              geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
#
# duppSorensen(allOncoGeneLists[["atlas"]], allOncoGeneLists[["sanger"]],
#              boot = TRUE,
#              listNames = c("atlas", "sanger"),
#              onto = "BP", GOLevel = 5,
#              geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
# etc.

# To save time, build the contingency table first, and then compute from it:
?buildEnrichTable
tab <- buildEnrichTable(allOncoGeneLists[["atlas"]], allOncoGeneLists[["sanger"]],
                        listNames = c("atlas", "sanger"),
                        onto = "BP", GOLevel = 5,
                        geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")

tab
equivTestSorensen(tab)
equivTestSorensen(tab, boot = TRUE)
dSorensen(tab)
seSorensen(tab)
duppSorensen(tab)
set.seed(123)
duppSorensen(tab, boot = TRUE)

# To perform from scratch all pairwise tests is even much slower:
allTests <- equivTestSorensen(allOncoGeneLists,
                              onto = "BP", GOLevel = 5,
                              geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
getPvalue(allTests, simplify = FALSE)
getPvalue(allTests)
p.adjust(getPvalue(allTests), method = "holm")

# All pairwise bootstrap tests from scratch (even more time consuming):
# set.seed(123)
# allBootTests <- equivTestSorensen(allOncoGeneLists,
#                                   boot = TRUE,
#                                   onto = "BP", GOLevel = 5,
#                                   geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
# Warnings are issued (not all bootstrap replicates may conduct to finite statistics):
# getNboot(allBootTests)

# It is much faster to upgrade the above tests to bootstrap tests:
set.seed(123)
allBootTests <- upgrade(allTests, boot = TRUE)
getPvalue(allBootTests, simplify = FALSE)

# To adjust for testing multiplicity:
# p.adjust(getPvalue(allBootTests), method = "holm")



# Obviously, building all pairwise contingency tables is even much slower
# allPairDiss <- dSorensen(allOncoGeneLists,
#                          onto = "BP", GOLevel = 5,
#                          geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
# allPairDiss

# Similarly, to obtain all pairwise...
# Standard errors:
# seSorensen(allOncoGeneLists,
#            onto = "BP", GOLevel = 5,
#            geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
#
# Upper confidence interval limits:
# duppSorensen(allOncoGeneLists,
#              onto = "BP", GOLevel = 5,
#              geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
# etc.

# Tipically, in a real study it would be interesting to scan tests
# along some ontologies and levels inside these ontologies:
# (which obviously will be a quite slow process)
# set.seed(123)
# allBootTests_BP_MF_lev4to8 <- allEquivTestSorensen(allOncoGeneLists,
#                                                    boot = TRUE,
#                                                    geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#                                                    ontos = c("BP", "MF"), GOLevels = 4:8)
# getPvalue(allBootTests_BP_MF_lev4to8)
# getNboot(allBootTests_BP_MF_lev4to8)
