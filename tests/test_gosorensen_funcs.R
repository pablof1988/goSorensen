library(goSorensen)

# A contingency table of GO terms mutual enrichment
# between gene lists "atlas" and "sanger":
data("cont_atlas.sanger_BP4")
cont_atlas.sanger_BP4
?cont_atlas.sanger_BP4
class(cont_atlas.sanger_BP4)

# Sorensen-Dice dissimilarity on this contingency table:
?dSorensen
dSorensen(cont_atlas.sanger_BP4)

# Standard error of this Sorensen-Dice dissimilarity estimate:
?seSorensen
seSorensen(cont_atlas.sanger_BP4)

# Upper 95% confidence limit for the Sorensen-Dice dissimilarity:
?duppSorensen
duppSorensen(cont_atlas.sanger_BP4)
# This confidence limit is based on an assimptotic normal N(0,1)
# approximation to the distribution of (dSampl - d) / se, where
# dSampl stands for the sample dissimilarity, d for the true dissimilarity
# and se for the sample dissimilarity standard error estimate.

# Upper confidence limit but using a Student's t instead of a N(0,1)
# (just as an example, not recommended -no theoretical justification)
df <- sum(cont_atlas.sanger_BP4[1:3]) - 2
duppSorensen(cont_atlas.sanger_BP4, z.conf.level = qt(1 - 0.95, df))

# Upper confidence limit but using a bootstrap approximation
# to the sampling distribution, instead of a N(0,1)
set.seed(123)
duppSorensen(cont_atlas.sanger_BP4, boot = TRUE)

# Some computations on diverse data structures:
badConti <- as.table(matrix(c(501, 27, 36, 12, 43, 15, 0, 0, 0),
                            nrow = 3, ncol = 3,
                            dimnames = list(c("a1","a2","a3"),
                                            c("b1", "b2","b3"))))
tryCatch(nice2x2Table(badConti), error = function(e) {return(e)})

incompleteConti <- badConti[1,1:min(2,ncol(badConti)), drop = FALSE]
incompleteConti
tryCatch(nice2x2Table(incompleteConti), error = function(e) {return(e)})

contiAsVector <- c(32, 21, 81, 1439)
nice2x2Table(contiAsVector)
contiAsVector.mat <- matrix(contiAsVector, nrow = 2)
contiAsVector.mat
contiAsVectorLen3 <- c(32, 21, 81)
nice2x2Table(contiAsVectorLen3)

tryCatch(dSorensen(badConti), error = function(e) {return(e)})

# Apparently, the next order works fine, but returns a wrong value!
dSorensen(badConti, check.table = FALSE)

tryCatch(dSorensen(incompleteConti), error = function(e) {return(e)})
dSorensen(contiAsVector)
dSorensen(contiAsVector.mat)
dSorensen(contiAsVectorLen3)
dSorensen(contiAsVectorLen3, check.table = FALSE)
dSorensen(c(0,0,0,45))

tryCatch(seSorensen(badConti), error = function(e) {return(e)})
tryCatch(seSorensen(incompleteConti), error = function(e) {return(e)})
seSorensen(contiAsVector)
seSorensen(contiAsVector.mat)
seSorensen(contiAsVectorLen3)
seSorensen(contiAsVectorLen3, check.table = FALSE)
tryCatch(seSorensen(contiAsVectorLen3, check.table = "not"), error = function(e) {return(e)})
seSorensen(c(0,0,0,45))

tryCatch(duppSorensen(badConti), error = function(e) {return(e)})
tryCatch(duppSorensen(incompleteConti), error = function(e) {return(e)})
duppSorensen(contiAsVector)
duppSorensen(contiAsVector.mat)
set.seed(123)
duppSorensen(contiAsVector, boot = TRUE)
set.seed(123)
duppSorensen(contiAsVector.mat, boot = TRUE)
duppSorensen(contiAsVectorLen3)
# Bootstrapping requires full contingency tables (4 values)
set.seed(123)
tryCatch(duppSorensen(contiAsVectorLen3, boot = TRUE), error = function(e) {return(e)})
duppSorensen(c(0,0,0,45))

# Equivalence test, H0: d >= d0 vs  H1: d < d0 (d0 = 0.4444)
?equivTestSorensen
equiv.atlas.sanger <- equivTestSorensen(cont_atlas.sanger_BP4)
equiv.atlas.sanger
getTable(equiv.atlas.sanger)
getPvalue(equiv.atlas.sanger)

tryCatch(equivTestSorensen(badConti), error = function(e) {return(e)})
tryCatch(equivTestSorensen(incompleteConti), error = function(e) {return(e)})
equivTestSorensen(contiAsVector)
equivTestSorensen(contiAsVector.mat)
set.seed(123)
equivTestSorensen(contiAsVector.mat, boot = TRUE)
equivTestSorensen(contiAsVectorLen3)

tryCatch(equivTestSorensen(contiAsVectorLen3, boot = TRUE), error = function(e) {return(e)})

equivTestSorensen(c(0,0,0,45))

# Sorensen-Dice computations from scratch, directly from gene lists
data(allOncoGeneLists)
?allOncoGeneLists

library(org.Hs.eg.db)
humanEntrezIDs <- keys(org.Hs.eg.db, keytype = "ENTREZID")
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
# To know the number of planned bootstrap replicates:
getNboot(bootTest)
# To know the number of valid bootstrap replicates:
getEffNboot(bootTest)

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
# set.seed(123)
# duppSorensen(allOncoGeneLists[["atlas"]], allOncoGeneLists[["sanger"]],
#              boot = TRUE,
#              listNames = c("atlas", "sanger"),
#              onto = "BP", GOLevel = 5,
#              geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
# etc.

# To build the contingency table first and then compute from it, may be a more flexible
# and saving time strategy, in general:
?buildEnrichTable
tab <- buildEnrichTable(allOncoGeneLists[["atlas"]], allOncoGeneLists[["sanger"]],
                        listNames = c("atlas", "sanger"),
                        onto = "BP", GOLevel = 5,
                        geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")

tab

# (Here, an obvious faster possibility would be to recover the enrichment contingency
# table from the previous normal test result:)
tab <- getTable(normTest)
tab

tst <- equivTestSorensen(tab)
tst
set.seed(123)
bootTst <- equivTestSorensen(tab, boot = TRUE)
bootTst

dSorensen(tab)
seSorensen(tab)
# or:
getDissimilarity(tst)

duppSorensen(tab)
getUpper(tst)

set.seed(123)
duppSorensen(tab, boot = TRUE)
getUpper(bootTst)

# To perform from scratch all pairwise tests (or other Sorensen-Dice computations)
# is even much slower. For example, all pairwise...
# Dissimilarities:
# # allPairDiss <- dSorensen(allOncoGeneLists,
# #                          onto = "BP", GOLevel = 5,
# #                          geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
# # allPairDiss
#
# # Still time consuming but faster: build all tables computing in parallel:
# allPairDiss <- dSorensen(allOncoGeneLists,
#                          onto = "BP", GOLevel = 5,
#                          geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#                          parallel = TRUE)
# allPairDiss

# Standard errors:
# seSorensen(allOncoGeneLists,
#            onto = "BP", GOLevel = 5,
#            geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
#
# Upper confidence interval limits:
# duppSorensen(allOncoGeneLists,
#              onto = "BP", GOLevel = 5,
#              geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
# All pairwise asymptotic normal tests:
# allTests <- equivTestSorensen(allOncoGeneLists,
#                               onto = "BP", GOLevel = 5,
#                               geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
# getPvalue(allTests, simplify = FALSE)
# getPvalue(allTests)
# p.adjust(getPvalue(allTests), method = "holm")
# To perform all pairwise bootstrap tests from scratch is (slightly)
# even more time consuming:
# set.seed(123)
# allBootTests <- equivTestSorensen(allOncoGeneLists,
#                                   boot = TRUE,
#                                   onto = "BP", GOLevel = 5,
#                                   geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
# Not all bootstrap replicates may conduct to finite statistics:
# getNboot(allBootTests)

# Given the normal tests (object 'allTests'), it is much faster to upgrade
# it to have the bootstrap tests:
# set.seed(123)
# allBootTests <- upgrade(allTests, boot = TRUE)
# getPvalue(allBootTests, simplify = FALSE)

# Again, the faster and more flexible possibility may be:
# 1) First, build all pairwise enrichment contingency tables (slow first step):
# allTabsBP.4 <- buildEnrichTable(allOncoGeneLists,
#                                 onto = "BP", GOLevel = 5,
#                                 geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
# allTabsBP.4

# Better, directly use the dataset available at this package, goSorensen:
data("cont_all_BP4")
cont_all_BP4
class(cont_all_BP4)
# 2) Then perform all required computatios from these enrichment contingency tables...
# All pairwise tests:
allTests <- equivTestSorensen(cont_all_BP4)
allTests
class(allTests)
set.seed(123)
allBootTests <- equivTestSorensen(cont_all_BP4, boot = TRUE)
allBootTests
class(allBootTests)
getPvalue(allBootTests, simplify = FALSE)
getEffNboot(allBootTests)

# To adjust for testing multiplicity:
p.adjust(getPvalue(allBootTests), method = "holm")

# If only partial statistics are desired:
dSorensen(cont_all_BP4)
duppSorensen(cont_all_BP4)
seSorensen(cont_all_BP4)


# Tipically, in a real study it would be interesting to scan tests
# along some ontologies and levels inside these ontologies:
# (which obviously will be a quite slow process)
# gc()
# set.seed(123)
# allBootTests_BP_MF_lev4to8 <- allEquivTestSorensen(allOncoGeneLists,
#                                                    boot = TRUE,
#                                                    geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#                                                    ontos = c("BP", "MF"), GOLevels = 4:8)
# getPvalue(allBootTests_BP_MF_lev4to8)
# getEffNboot(allBootTests_BP_MF_lev4to8)
