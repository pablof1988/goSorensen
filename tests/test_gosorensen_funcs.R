library(goSorensen)

# A contingency table of GO terms mutual enrichment
# between gene lists "atlas" and "sanger":
data(tab_atlas.sanger_BP3)
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

# Apparently this works fine, but returns a wrong value:
dSorensen(badConti, check.table = FALSE)

tryCatch(dSorensen(incompleteConti), error = function(e) {return(e)})
dSorensen(contiAsVector)
dSorensen(contiAsVector.mat)
dSorensen(contiAsVectorLen3)
dSorensen(contiAsVectorLen3, check.table = FALSE)

tryCatch(seSorensen(badConti), error = function(e) {return(e)})
tryCatch(seSorensen(incompleteConti), error = function(e) {return(e)})
seSorensen(contiAsVector)
seSorensen(contiAsVector.mat)
seSorensen(contiAsVectorLen3)
seSorensen(contiAsVectorLen3, check.table = FALSE)

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

# Equivalence test, H0: d >= d0 vs  H1: d < d0 (d0 = 0.4444)
?equivTestSorensen
equiv.atlas.sanger <- equivTestSorensen(tab_atlas.sanger_BP3)
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

# Sorensen-Dice computations from scratch, directly from gene lists
data(allOncoGeneLists)
?allOncoGeneLists
# Gene universe:
data(humanEntrezIDs)
# First, the mutual GO node enrichment tables are built, then computations
# proceed from these contingency tables.
# Building the contingency tables is a slow process (many enrichment tests)
normTest <- equivTestSorensen(allOncoGeneLists[["atlas"]], 
                              allOncoGeneLists[["sanger"]],
                              listNames = c("atlas", "sanger"),
                              onto = "BP", GOLevel = 5,
                              geneUniverse = humanEntrezIDs, 
                              orgPackg = "org.Hs.eg.db")
normTest

# It is much faster to upgrade 'normTest' to be a bootstrap test:
set.seed(123)
bootTest <- upgrade(normTest, boot = TRUE)
bootTest
# To know the number of valid bootstrap replicates:
getNboot(bootTest)


# To save time, build the contingency table first, and then compute from it:
?buildEnrichTable
tab <- buildEnrichTable(allOncoGeneLists[["atlas"]], 
                        allOncoGeneLists[["sanger"]],
                        listNames = c("atlas", "sanger"),
                        onto = "BP", GOLevel = 5,
                        geneUniverse = humanEntrezIDs, 
                        orgPackg = "org.Hs.eg.db")

tab
equivTestSorensen(tab)
equivTestSorensen(tab, boot = TRUE)
dSorensen(tab)
seSorensen(tab)
duppSorensen(tab)
set.seed(123)
duppSorensen(tab, boot = TRUE)