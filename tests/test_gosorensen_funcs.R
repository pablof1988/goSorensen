library(goSorensen)
library(equivStandardTest)

data(humanEntrezIDs)

# ------ontology BP level 3 ----------------- atlas , sanger ------------------
#                    Enriched in sanger
# Enriched in atlas  FALSE TRUE
#              FALSE   471    1
#              TRUE     30   56
tab_atlas.sanger_BP3
?tab_atlas.sanger_BP3
class(tab_atlas.sanger_BP3)

dSorensen(tab_atlas.sanger_BP3)
seSorensen(tab_atlas.sanger_BP3)
duppSorensen(tab_atlas.sanger_BP3)
equivTestSorensen(tab_atlas.sanger_BP3)


badConti <- as.table(matrix(c(501, 27, 36, 12, 43, 15, 0, 0, 0),
                            nrow = 3, ncol = 3,
                            dimnames = list(c("a1","a2","a3"),
                                            c("b1", "b2","b3"))))
nice2x2Table(badConti)
nice2x2Table(badConti, listNames = c("a gene list","another"))
incompleteConti <- badConti[1,1:min(2,ncol(badConti)), drop = FALSE]
incompleteConti
nice2x2Table(incompleteConti)
nice2x2Table(incompleteConti, listNames = c("a gene list","another"))

anotherIncomplete <- matrix(c(210, 12), ncol = 1, nrow = 2)
anotherIncomplete
nice2x2Table(anotherIncomplete)

contiAsVector <- c(1439, 32, 21, 81)
nice2x2Table(contiAsVector)
contiAsVector.mat <- matrix(contiAsVector, nrow = 2)
contiAsVector.mat
contiAsVectorLen3 <- c(32, 21, 81)
nice2x2Table(contiAsVectorLen3)
nice2x2Table(contiAsVectorLen3, n = 1573)
nice2x2Table(contiAsVectorLen3, n = 1000)
try(nice2x2Table(contiAsVectorLen3, n = 10))

?tab_atlas.sanger_BP3
dSorensen(tab_atlas.sanger_BP3)

dSorensen(badConti)
dSorensen(badConti, check.table = FALSE)
# Wrong value!
# Better:
dSorensen(nice2x2Table(badConti), check.table = FALSE)

dSorensen(incompleteConti)
dSorensen(anotherIncomplete)
dSorensen(contiAsVector)
dSorensen(contiAsVector.mat)
dSorensen(contiAsVectorLen3)
dSorensen(contiAsVectorLen3, check.table = FALSE)

?pbtGeneLists
dSorensen(pbtGeneLists[[2]], pbtGeneLists[[4]],
          listNames = names(pbtGeneLists)[c(2,4)],
          onto = "BP", GOLevel = 5,
          geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
# Essentially, the above code makes the same as:
pbtBP5.IRITD3vsKT1 <- nice2x2Table(
  crossTabGOIDs4GeneLists(pbtGeneLists[[2]], pbtGeneLists[[4]],
                          onto = "BP", GOLevel = 5,
                          geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db"),
  listNames = names(pbtGeneLists)[c(2,4)]
)
dSorensen(pbtBP5.IRITD3vsKT1)

dSorensen(pbtGeneLists,
          onto = "BP", GOLevel = 5,
          geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")

seSorensen(badConti)
seSorensen(incompleteConti)
seSorensen(anotherIncomplete)
seSorensen(contiAsVector)
seSorensen(contiAsVector.mat)
seSorensen(contiAsVectorLen3)
seSorensen(contiAsVectorLen3, check.table = FALSE)
# To compute 'se' for a proportions table:
seSorensen(contiAsVector.mat/sum(contiAsVector.mat), n = sum(contiAsVector.mat[2:4]))
# # Alternatively:
# seSorensen(contiAsVector.mat/sum(contiAsVector.mat)) / sqrt(sum(contiAsVector.mat))
seSorensen(contiAsVectorLen3 / sum(contiAsVectorLen3), n = sum(contiAsVectorLen3))

seSorensen(pbtGeneLists[[2]], pbtGeneLists[[4]],
           listNames = names(pbtGeneLists)[c(2,4)],
           onto = "BP", GOLevel = 5,
           geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
# # Essentially, the above code makes the same as:
# pbtBP5.IRITD3vsKT1 <- nice2x2Table(
#   crossTabGOIDs4GeneLists(pbtGeneLists[[2]], pbtGeneLists[[4]],
#                           onto = "BP", GOLevel = 5,
#                           geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db"),
#   listNames = names(pbtGeneLists)[c(2,4)]
# )
# seSorensen(pbtBP5.IRITD3vsKT1)

seSorensen(pbtGeneLists,
           onto = "BP", GOLevel = 5,
           geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")


duppSorensen(badConti)
duppSorensen(incompleteConti)
duppSorensen(anotherIncomplete)
duppSorensen(contiAsVector)
duppSorensen(contiAsVector.mat)
duppSorensen(contiAsVectorLen3)
duppSorensen(contiAsVectorLen3, check.table = FALSE)
duppSorensen(contiAsVectorLen3 / sum(contiAsVectorLen3), n = sum(contiAsVectorLen3))

duppSorensen(pbtGeneLists[[2]], pbtGeneLists[[4]],
             listNames = names(pbtGeneLists)[c(2,4)],
             onto = "BP", GOLevel = 5,
             geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
# # Essentially, the above code makes the same as:
# pbtBP5.IRITD3vsKT1 <- nice2x2Table(
#   crossTabGOIDs4GeneLists(pbtGeneLists[[2]], pbtGeneLists[[4]],
#                           onto = "BP", GOLevel = 5,
#                           geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db"),
#   listNames = names(pbtGeneLists)[c(2,4)]
# )
# duppSorensen(pbtBP5.IRITD3vsKT1)

duppSorensen(pbtGeneLists,
             onto = "BP", GOLevel = 5,
             geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")


equivTestSorensen(tab_atlas.sanger_BP3)
equivTestSorensen(badConti)
equivTestSorensen(incompleteConti)
equivTestSorensen(anotherIncomplete)
equivTestSorensen(contiAsVector)
equivTestSorensen(contiAsVector.mat)
equivTestSorensen(contiAsVector.mat / sum(contiAsVector.mat), n = sum(contiAsVector.mat[2:4]))
equivTestSorensen(contiAsVectorLen3)
equivTestSorensen(contiAsVectorLen3, check.table = FALSE)
equivTestSorensen(contiAsVectorLen3 / sum(contiAsVectorLen3), n = sum(contiAsVectorLen3))

equivTestSorensen(pbtGeneLists[[2]], pbtGeneLists[[4]],
                  listNames = names(pbtGeneLists)[c(2,4)],
                  onto = "BP", GOLevel = 5,
                  geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
# # Essentially, the above code makes the same as:
# pbtBP5.IRITD3vsKT1 <- nice2x2Table(
#   crossTabGOIDs4GeneLists(pbtGeneLists[[2]], pbtGeneLists[[4]],
#                           onto = "BP", GOLevel = 5,
#                           geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db"),
#   listNames = names(pbtGeneLists)[c(2,4)]
# )
# equivTestSorensen(pbtBP5.IRITD3vsKT1)

equivTestSorensen(pbtGeneLists,
                  onto = "BP", GOLevel = 5,
                  geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
