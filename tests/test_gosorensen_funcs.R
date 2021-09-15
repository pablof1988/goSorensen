library(goSorensen)
library(equivStandardTest)

data(humanEntrezIDs)

tab_atlas.sanger_BP3
?tab_atlas.sanger_BP3
class(tab_atlas.sanger_BP3)

?dSorensen
dSorensen(tab_atlas.sanger_BP3)
?seSorensen
seSorensen(tab_atlas.sanger_BP3)
?duppSorensen
duppSorensen(tab_atlas.sanger_BP3)
# Upper confidence limit but using a Student's t instead of a N(0,1)
df <- sum(tab_atlas.sanger_BP3[1:3]) - 2
duppSorensen(tab_atlas.sanger_BP3, z.conf.level = qt(1 - 0.95, df))

qBoot <- function(p, xTab, nboot = 10000) {
  n <- sum(xTab)
  pTab <- xTab / n
  # pEnrich <- sum(pTab[1:3])
  # pTab <- pTab[1:3] / pEnrich
  dS <- dSorensen(xTab)
  bootTabs <- rmultinom(nboot, size = n, prob = pTab)
  # nEnrich <- rbinom(nboot, n, pEnrich)
  tStats <- apply(bootTabs, 2, function(xBoot) {
    xBoot <- matrix(xBoot, nrow = 2)
    (dSorensen(xBoot) - dS) / seSorensen(xBoot)
  })
  return(quantile(tStats, probs = p))
}

set.seed(123)
qBoot(1 - 0.95, tab_atlas.sanger_BP3)

# Upper confidence limit but using a bootstrap approximation
# to the sampling distribution, instead of a N(0,1)
set.seed(123)
duppSorensen(tab_atlas.sanger_BP3, z.conf.level = qBoot(1 - 0.95, tab_atlas.sanger_BP3))
# The same but with an internal bootstrap procedure:
set.seed(123)
duppSorensen(tab_atlas.sanger_BP3, boot = TRUE)

?equivTestSorensen
testResult <- equivTestSorensen(tab_atlas.sanger_BP3)
getPvalue(testResult)
getTable(testResult)

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

dSorensen(badConti)
dSorensen(badConti, check.table = FALSE)
# Wrong value!

dSorensen(incompleteConti)
dSorensen(contiAsVector)
dSorensen(contiAsVector.mat)
dSorensen(contiAsVectorLen3)
dSorensen(contiAsVectorLen3, check.table = FALSE)
# dSorensen(32)
# dSorensen(32, n = 15)
# dSorensen(32, n = 134)
# dSorensen(32/134)

?pbtGeneLists
dSorensen(pbtGeneLists[[2]], pbtGeneLists[[4]],
          listNames = names(pbtGeneLists)[c(2,4)],
          onto = "BP", GOLevel = 5,
          geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")

dSorensen(pbtGeneLists,
          onto = "BP", GOLevel = 5,
          geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")

seSorensen(badConti)
seSorensen(incompleteConti)
seSorensen(contiAsVector)
seSorensen(contiAsVector.mat)
seSorensen(contiAsVectorLen3)
seSorensen(contiAsVectorLen3, check.table = FALSE)
# # To compute 'se' for a proportions table:
# seSorensen(contiAsVector.mat/sum(contiAsVector.mat[1:3]), n = sum(contiAsVector.mat[1:3]))
# seSorensen(contiAsVectorLen3 / sum(contiAsVectorLen3), n = sum(contiAsVectorLen3))

seSorensen(pbtGeneLists[[2]], pbtGeneLists[[4]],
           listNames = names(pbtGeneLists)[c(2,4)],
           onto = "BP", GOLevel = 5,
           geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")

seSorensen(pbtGeneLists,
           onto = "BP", GOLevel = 5,
           geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")


duppSorensen(badConti)
duppSorensen(incompleteConti)
duppSorensen(contiAsVector)
duppSorensen(contiAsVector.mat)
set.seed(123)
duppSorensen(contiAsVector, boot = TRUE)
set.seed(123)
duppSorensen(contiAsVector.mat, boot = TRUE)

duppSorensen(contiAsVectorLen3)
set.seed(123)
duppSorensen(contiAsVectorLen3, boot = TRUE)
duppSorensen(contiAsVectorLen3, check.table = FALSE)
# duppSorensen(contiAsVectorLen3 / sum(contiAsVectorLen3), n = sum(contiAsVectorLen3))

duppSorensen(pbtGeneLists[[2]], pbtGeneLists[[4]],
             listNames = names(pbtGeneLists)[c(2,4)],
             onto = "BP", GOLevel = 5,
             geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")

duppSorensen(pbtGeneLists,
             onto = "BP", GOLevel = 5,
             geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")

duppSorensen(pbtGeneLists[[2]], pbtGeneLists[[4]],
             boot = TRUE,
             listNames = names(pbtGeneLists)[c(2,4)],
             onto = "BP", GOLevel = 5,
             geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")

duppSorensen(pbtGeneLists,
             boot = TRUE,
             onto = "BP", GOLevel = 5,
             geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")

duppSorensen(pbtGeneLists[["KT1.1"]], pbtGeneLists[["ENDAT"]],
             boot = TRUE,
             listNames = names(pbtGeneLists)[c(2,4)],
             onto = "BP", GOLevel = 5,
             geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")

tab <- buildEnrichTable(pbtGeneLists[["KT1.1"]], pbtGeneLists[["ENDAT"]],
                 listNames = names(pbtGeneLists)[c(2,4)],
                 onto = "BP", GOLevel = 5,
                 geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")


equiv.atlas.sanger <- equivTestSorensen(tab_atlas.sanger_BP3)
equiv.atlas.sanger
getTable(equiv.atlas.sanger)
getPvalue(equiv.atlas.sanger)

equivTestSorensen(badConti)
equivTestSorensen(incompleteConti)
equivTestSorensen(contiAsVector)
equivTestSorensen(contiAsVector.mat)
equivTestSorensen(contiAsVector.mat, boot = TRUE)
equivTestSorensen(contiAsVectorLen3)
equivTestSorensen(contiAsVectorLen3, boot = TRUE)

equivTestSorensen(allOncoGeneLists[[2]], allOncoGeneLists[[4]],
                  listNames = names(allOncoGeneLists)[c(2,4)],
                  onto = "BP", GOLevel = 5,
                  geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")

set.seed(123)
equivTestSorensen(allOncoGeneLists[[2]], allOncoGeneLists[[4]],
                  listNames = names(allOncoGeneLists)[c(2,4)],
                  boot = TRUE,
                  onto = "BP", GOLevel = 5,
                  geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")

allTests <- equivTestSorensen(allOncoGeneLists,
                              onto = "BP", GOLevel = 5,
                              geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
getPvalue(allTests)
p.adjust(getPvalue(allTests), method = "holm")

equivTestSorensen(allOncoGeneLists,
                  boot = TRUE,
                  onto = "BP", GOLevel = 5,
                  geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")

allEquivTestSorensen(allOncoGeneLists,
                     boot = TRUE,
                     geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
                     ontos = c("BP", "MF"), GOLevels = 4:5)
