% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getTable.R
\name{getTable}
\alias{getTable}
\alias{getTable.equivSDhtest}
\alias{getTable.equivSDhtestList}
\alias{getTable.AllEquivSDhtest}
\title{Access to the contingency table of mutual enrichment of one or more equivalence test results}
\usage{
getTable(x, ...)

\method{getTable}{equivSDhtest}(x, ...)

\method{getTable}{equivSDhtestList}(x, ...)

\method{getTable}{AllEquivSDhtest}(x, onto, GOLevel, listNames, ...)
}
\arguments{
\item{x}{an object of class "equivSDhtest" or "equivSDhtestList" or "allEquivSDtest".}

\item{...}{Additional parameters.}

\item{onto}{character, a vector with one or more of "BP", "CC" or "MF", ontologies to access.}

\item{GOLevel}{numeric or character, a vector with one or more GO levels to access.
See the details section and the examples.}

\item{listNames}{character(2), the names of a pair of gene lists.}
}
\value{
An object of class "table", the 2x2 enrichment contingeny table of mutual enrichment
in two gene lists, built to perform the equivalence test based on the Sorensen-Dice dissimilarity.

When \code{x} is an object of class "equivSDhtest" (i.e., the result of a single
equivalence test), the returned value is an object of class "table", the 2x2 enrichment
contingeny table of mutual enrichment in two gene lists, built to perform the equivalence test
based on the Sorensen-Dice dissimilarity.
For an object of class "equivSDhtestList" (i.e. all pairwise tests for a
set of gene lists), the resulting value is a list with all the tables built in all those
tests. If \code{x} is an object of class "allEquivSDtest"
(i.e., the test iterated along GO ontologies and levels), the preceding result is returned
as a list along the ontologies, levels and pairs of gene lists specified by the arguments
\code{onto, GOlevel} and \code{listNames} (or all ontologies, levels or pairs of gene lists
present in \code{x} if one or more of these arguments are missing).
}
\description{
Given objects representing the result(s) of one or more equivalence tests
(classes "equivSDhtest", "equivSDhtestList" or "allEquivSDtest", i.e., the
result of functions 'equivTestSorensen' and 'allEquivTestSorensen')
this function returns the contingency tables from which the tests were performed.
}
\details{
Argument \code{GOLevel} can be of class "character" or "numeric". In the first case, the GO
levels must be specified like \code{"level 6"} or \code{c("level 4", "level 5", "level 6")}
In the second case ("numeric"), the GO levels must be specified like\code{6} or \code{4:6}.
}
\section{Methods (by class)}{
\itemize{
\item \code{getTable(equivSDhtest)}: S3 method for class "equivSDhtest"

\item \code{getTable(equivSDhtestList)}: S3 method for class "equivSDhtestList"

\item \code{getTable(AllEquivSDhtest)}: S3 method for class "AllEquivSDhtest"

}}
\examples{
# Dataset 'allOncoGeneLists' contains the result of the equivalence test between gene lists
# 'sanger' and 'atlas', at level 4 of the BP ontology:
data(eqTest_atlas.sanger_BP4)
eqTest_atlas.sanger_BP4
class(eqTest_atlas.sanger_BP4)
# This may correspond to the result of code like:
# eqTest_atlas.sanger_BP4 <- equivTestSorensen(
#   allOncoGeneLists[["sanger"]], allOncoGeneLists[["atlas"]],
#   geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#   onto = "BP", GOLevel = 4, listNames = c("sanger", "atlas"))
# (But results may vary according to GO updating)
getTable(eqTest_atlas.sanger_BP4)

# All pairwise equivalence tests at level 4 of the BP ontology
data(eqTest_all_BP4)
?eqTest_all_BP4
class(eqTest_all_BP4)
# This may correspond to a call like:
# eqTest_all_BP4 <- equivTestSorensen(allOncoGeneLists,
#                           geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#                           onto = "BP", GOLevel = 4)
getTable(eqTest_all_BP4)

# Equivalence test iterated over all GO ontologies and levels 3 to 10:
data(allEqTests)
?allEqTests
class(allEqTests)
# This may correspond to code like:
# allEqTests <- allEquivTestSorensen(allOncoGeneLists,
#                                             geneUniverse = humanEntrezIDs,
#                                             orgPackg = "org.Hs.eg.db")
# (By default, the tests are iterated over all GO ontologies and for levels 3 to 10)
# All 2x2 contingecy tables of joint enrichment:
getTable(allEqTests)
# Contingency tables only for some GO ontologies, levels or pairs of gene lists:
getTable(allEqTests, GOLevel = "level 6")
getTable(allEqTests, GOLevel = 6)
getTable(allEqTests, GOLevel = seq.int(4,6), listNames = c("atlas", "sanger"))
getTable(allEqTests, GOLevel = "level 6", onto = "BP")
getTable(allEqTests, GOLevel = "level 6", onto = "BP",
         listNames = c("waldman", "sanger"))


}
