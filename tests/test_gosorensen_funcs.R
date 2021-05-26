library(goSorensen)
library(equivStandardTest)

data(humanEntrezIDs)

# ------ontology BP level 3 ----------------- atlas , sanger ------------------
#   Enriched in sanger
# Enriched in atlas FALSE TRUE
# FALSE   471    1
# TRUE     30   56
tab_atlas.sanger_BP3
?tab_atlas.sanger_BP3
class(tab_atlas.sanger_BP3)

dSorensen(tab_atlas.sanger_BP3)
seSorensen(tab_atlas.sanger_BP3)
duppSorensen(tab_atlas.sanger_BP3)
equivTestSorensen(tab_atlas.sanger_BP3)


conti <- as.table(matrix(c(501, 27, 36, 12, 43, 15, 0, 0, 0),
                         nrow = 3, ncol = 3,
                         dimnames = list(c("a1","a2","a3"),
                                         c("b1", "b2","b3"))))
nice2x2Table(conti)
nice2x2Table(conti, listNames = c("a gene list","another"))
conti2 <- conti[1,1:min(2,ncol(conti)), drop = FALSE]
conti2
nice2x2Table(conti2)
nice2x2Table(conti2, listNames = c("a gene list","another"))

conti3 <- matrix(c(210, 12), ncol = 2, nrow = 1)
conti3
nice2x2Table(conti3)

conti4 <- c(1439, 32, 21, 81)
nice2x2Table(conti4)
conti4.mat <- matrix(conti4, nrow = 2)
conti4.mat
conti5 <- c(32, 21, 81)
nice2x2Table(conti5, n = 1573)
nice2x2Table(conti5, n = 1000)
try(nice2x2Table(conti5, n = 10))

?tab_atlas.sanger_BP3
dSorensen(tab_atlas.sanger_BP3)

dSorensen(conti)
dSorensen(conti, check.table = FALSE)
# Wrong value!
# Better:
dSorensen(nice2x2Table(conti), check.table = FALSE)

dSorensen(conti2)
dSorensen(conti3)
dSorensen(conti4)
dSorensen(conti4.mat)
dSorensen(conti5, n = 1573)
dSorensen(conti5, n = 1000)
dSorensen(conti5, n = 1573, check.table = FALSE)

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

seSorensen(conti)
seSorensen(conti2)
seSorensen(conti3)
seSorensen(conti4)
seSorensen(conti4.mat)
seSorensen(conti5, n = 1573)
seSorensen(conti5, n = 1000)
seSorensen(conti5, n = 1573, check.table = FALSE)
# To compute 'se' for a proportions table:
seSorensen(conti4.mat/sum(conti4.mat)) / sqrt(sum(conti4.mat))
# In general:
# seSorensen(pij) / sqrt(n)

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
seSorensen(pbtBP5.IRITD3vsKT1)

seSorensen(pbtGeneLists,
           onto = "BP", GOLevel = 5,
           geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")


duppSorensen(conti)
duppSorensen(conti2)
duppSorensen(conti3)
duppSorensen(conti4)
duppSorensen(conti4.mat)
duppSorensen(conti5, n = 1573)
duppSorensen(conti5, n = 1000)
duppSorensen(conti5, n = 1573, check.table = FALSE)

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
duppSorensen(pbtBP5.IRITD3vsKT1)

duppSorensen(pbtGeneLists,
             onto = "BP", GOLevel = 5,
             geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")


equivTestSorensen(tab_atlas.sanger_BP3)
equivTestSorensen(conti)
equivTestSorensen(conti2)
equivTestSorensen(conti3)
equivTestSorensen(conti4)
equivTestSorensen(conti4.mat)
equivTestSorensen(conti5, n = 1573)
equivTestSorensen(conti5, n = 1000)
equivTestSorensen(conti5, n = 1573, check.table = FALSE)

# To compute se for a proportions table:
# seSorensen(prova.tab/sum(prova.tab)) / sqrt(sum(prova.tab))
# In general:
# seSorensen(pij) / sqrt(n)

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
equivTestSorensen(pbtBP5.IRITD3vsKT1)

equivTestSorensen(pbtGeneLists,
                  onto = "BP", GOLevel = 5,
                  geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
