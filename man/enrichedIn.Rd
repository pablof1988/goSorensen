% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrichedIn.R
\name{enrichedIn}
\alias{enrichedIn}
\alias{enrichedIn.default}
\alias{enrichedIn.character}
\alias{enrichedIn.list}
\title{This function builds a cross-tabulation of enriched (TRUE) and non-enriched (FALSE) GO terms vs. gene lists}
\usage{
enrichedIn(x, ...)

\method{enrichedIn}{default}(
  x,
  geneUniverse,
  orgPackg,
  onto,
  GOLevel,
  pAdjustMeth = "BH",
  pvalCutoff = 0.01,
  qvalCutoff = 0.05,
  parallel = FALSE,
  nOfCores = 1,
  onlyEnriched = TRUE,
  ...
)

\method{enrichedIn}{character}(
  x,
  geneUniverse,
  orgPackg,
  onto,
  GOLevel,
  pAdjustMeth = "BH",
  pvalCutoff = 0.01,
  qvalCutoff = 0.05,
  parallel = FALSE,
  nOfCores = 1,
  onlyEnriched = TRUE,
  ...
)

\method{enrichedIn}{list}(
  x,
  geneUniverse,
  orgPackg,
  onto,
  GOLevel,
  pAdjustMeth = "BH",
  pvalCutoff = 0.01,
  qvalCutoff = 0.05,
  parallel = FALSE,
  nOfCores = min(detectCores() - 1, length(x)),
  onlyEnriched = TRUE,
  ...
)
}
\arguments{
\item{x}{either an object of class "character" (or coerzable to "character") or "list".
In the "character" interface, these values should represent Entrez gene (or, in general, feature)
identifiers.
In the "list" interface, each element of the list must be a "character" vector of
Entrez identifiers}

\item{...}{Additional parameters}

\item{geneUniverse}{character vector containing the universe of genes from where gene lists have been extracted. This vector must be obtained from the annotation package declared in \code{orgPackg}. For more details, refer to vignette \href{../doc/goSorensen_Introduction.html}{goSorensen_Introduction}.}

\item{orgPackg}{A string with the name of the genomic annotation package corresponding to a specific species to be analyzed, which must be previously installed and activated. For more details, refer to vignette \href{../doc/goSorensen_Introduction.html}{goSorensen_Introduction}.}

\item{onto}{string describing the ontology. Belongs to c('BP', 'MF', 'CC')}

\item{GOLevel}{GO level, an integer}

\item{pAdjustMeth}{string describing the adjust method. Belongs to c('BH', 'BY', 'Bonf')}

\item{pvalCutoff}{adjusted pvalue cutoff on enrichment tests to report}

\item{qvalCutoff}{qvalue cutoff on enrichment tests to report as significant.
Tests must pass i) pvalueCutoff on unadjusted pvalues, ii) pvalueCutoff on adjusted pvalues and
iii) qvalueCutoff on qvalues to be reported}

\item{parallel}{Logical. Only in "list" interface. Defaults to FALSE, put it at TRUE for parallel computation}

\item{nOfCores}{Number of cores for parallel computations. Only in "list" interface}

\item{onlyEnriched}{logical. If TRUE (the default), the returned result only contains those GO terms
enriched in almost one of the gene lists}
}
\value{
In the "character" interface, a length k vector of TRUE/FALSE values
corresponding to enrichment or not of the GO terms at level 'GOLev' in ontology 'onto'.
If 'onlyEnriched' is FALSE, k corresponds to the total number of these GO terms. If 'onlyEnriched'
is TRUE (default) k is the number of enriched GO terms (and then all values in the resulting
vector are TRUE).
In the "list" interface, a logical matrix of TRUE/FALSE values indicating enrichment or not,
with k rows and s columns. s is the number of gene lists (the length of "list" 'x').
If 'onlyEnriched' is FALSE, k corresponds to the total number of GO terms at level 'GOLev' in ontology
'onto'. If 'onlyEnriched' is TRUE (default), the resulting matrix only contains the k rows
corresponding to GO terms enriched in almost one of these s gene lists.
In both interfaces ("character" or "list"), the result also has an attribute (\code{nTerms}) with
the total number of GO terms at level 'GOLev' in ontology 'onto'.
}
\description{
This function builds a cross-tabulation of enriched (TRUE) and non-enriched (FALSE) GO terms vs. gene lists
}
\details{
When the function argument \code{onlyEnriched} is FALSE, commonly the result is a sparse
but very large object. This function is primarily designed for internal use of function
\code{buildEnrichTable}, with argument \code{onlyEnriched} always put at its default TRUE value.
Then calls to \code{enrichedIn} result in much more compact objects, in general.

Argument \code{parallel} only applies to interface "list". Its default value ís "FALSE" and you
may consider the trade of between the time spent in initializing parallelization and the possible
time gain when parallelizing. It is difficult to establish a general guideline, but parallelizing
is only worthwhile when analyzing many gene lists, on the order of 30 or more, although it depends
a lot on each processor.
}
\section{Methods (by class)}{
\itemize{
\item \code{enrichedIn(default)}: S3 default method

\item \code{enrichedIn(character)}: S3 method for class "character"

\item \code{enrichedIn(list)}: S3 method for class "list"

}}
\examples{
# Obtaining ENTREZ identifiers for the gene universe of humans:
library(org.Hs.eg.db)
humanEntrezIDs <- keys(org.Hs.eg.db, keytype = "ENTREZID")

# Gene lists to be explored for enrichment:
data(allOncoGeneLists)
?allOncoGeneLists

# Computing the cross table:
enrichd <- enrichedIn(allOncoGeneLists[["Vogelstein"]],
                      geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
                      onto = "MF", GOLevel = 6)
enrichd

# Cross table of enriched GO terms (GO ontology MF, level 6) for all gene
# lists in 'allOncoGeneLists':
enrichedAllOncoMF.6 <- enrichedIn(allOncoGeneLists,
                          geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
                          onto = "MF", GOLevel = 6)
enrichedAllOncoMF.6
object.size(enrichedAllOncoMF.6)
# How many GO terms were tested for enrichment at ontology MF and level 6:
attr(enrichedAllOncoMF.6, "nTerms")
}
