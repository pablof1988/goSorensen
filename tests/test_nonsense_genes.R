library(goSorensen)

testError <- function(e) {return(e)}

tryCatch(dSorensen("Sec1", onto = "BP"), error = testError)

data(allOncoGeneLists)
?allOncoGeneLists
data(humanEntrezIDs)

# Non-sense random gene lists. Generating Entrez-like gene identifiers, but random:
set.seed(1234567)
genList1 <- unique(as.character(sample.int(99999, size = 100)))
genList2 <- unique(as.character(sample.int(99999, size = 100)))
# Gene identifiers are numbers like Entrez identifiers at 'humanEntrezIDs', but random.
dSorensen(genList1, genList2,
          listNames = c("genList1", "genList2"),
          onto = "BP", GOLevel = 5,
          geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
duppSorensen(genList1, genList2,
             listNames = c("genList1", "genList2"),
             onto = "BP", GOLevel = 5,
             geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
seSorensen(genList1, genList2,
           listNames = c("genList1", "genList2"),
           onto = "BP", GOLevel = 5,
           geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
nonSenseTst <- equivTestSorensen(genList1, genList2,
                                 listNames = c("genList1", "genList2"),
                                 onto = "BP", GOLevel = 5,
                                 geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
nonSenseTst
tab <- getTable(nonSenseTst)
tab
# Or, alternatively:
tab <- buildEnrichTable(genList1, genList2,
                        listNames = c("genList1", "genList2"),
                        onto = "BP", GOLevel = 5,
                        geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
tab
dSorensen(tab)
duppSorensen(tab)
seSorensen(tab)
equivTestSorensen(tab)

# Even more non-sense, letters non numeric-style like those at 'humanEntrezIDs':
set.seed(1234567)
genList1 <- unique(vapply(seq_len(100), function(i) {
  paste0(sample(c(letters, LETTERS), 6, replace = TRUE), collapse = "")
}, FUN.VALUE = character(1)))
genList2 <- unique(vapply(seq_len(100), function(i) {
  paste0(sample(c(letters, LETTERS), 6, replace = TRUE), collapse = "")
}, FUN.VALUE = character(1)))

# Gene identifiers incompatible with those at 'humanEntrezIDs':
dSorensen(genList1, genList2,
          listNames = c("genList1", "genList2"),
          onto = "BP", GOLevel = 5,
          geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
duppSorensen(genList1, genList2,
             listNames = c("genList1", "genList2"),
             onto = "BP", GOLevel = 5,
             geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
seSorensen(genList1, genList2,
           listNames = c("genList1", "genList2"),
           onto = "BP", GOLevel = 5,
           geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
nonSenseTst <- equivTestSorensen(genList1, genList2,
                                 listNames = c("genList1", "genList2"),
                                 onto = "BP", GOLevel = 5,
                                 geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
nonSenseTst
tab <- getTable(nonSenseTst)
tab
# Or, alternatively:
tab <- buildEnrichTable(genList1, genList2,
                        listNames = c("genList1", "genList2"),
                        onto = "BP", GOLevel = 5,
                        geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
tab
dSorensen(tab)
duppSorensen(tab)
seSorensen(tab)
equivTestSorensen(tab)
