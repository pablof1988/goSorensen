% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/duppsorensen.R
\name{duppSorensen}
\alias{duppSorensen}
\alias{duppSorensen.table}
\alias{duppSorensen.matrix}
\alias{duppSorensen.numeric}
\alias{duppSorensen.character}
\alias{duppSorensen.list}
\alias{duppSorensen.tableList}
\title{Upper limit of a one-sided confidence interval (0, dUpp] for the Sorensen-Dice dissimilarity}
\usage{
duppSorensen(x, ...)

\method{duppSorensen}{table}(
  x,
  dis = dSorensen.table(x, check.table = FALSE),
  se = seSorensen.table(x, check.table = FALSE),
  conf.level = 0.95,
  z.conf.level = qnorm(1 - conf.level),
  boot = FALSE,
  nboot = 10000,
  check.table = TRUE,
  ...
)

\method{duppSorensen}{matrix}(
  x,
  dis = dSorensen.matrix(x, check.table = FALSE),
  se = seSorensen.matrix(x, check.table = FALSE),
  conf.level = 0.95,
  z.conf.level = qnorm(1 - conf.level),
  boot = FALSE,
  nboot = 10000,
  check.table = TRUE,
  ...
)

\method{duppSorensen}{numeric}(
  x,
  dis = dSorensen.numeric(x, check.table = FALSE),
  se = seSorensen.numeric(x, check.table = FALSE),
  conf.level = 0.95,
  z.conf.level = qnorm(1 - conf.level),
  boot = FALSE,
  nboot = 10000,
  check.table = TRUE,
  ...
)

\method{duppSorensen}{character}(
  x,
  y,
  conf.level = 0.95,
  boot = FALSE,
  nboot = 10000,
  check.table = TRUE,
  ...
)

\method{duppSorensen}{list}(
  x,
  conf.level = 0.95,
  boot = FALSE,
  nboot = 10000,
  check.table = TRUE,
  ...
)

\method{duppSorensen}{tableList}(
  x,
  conf.level = 0.95,
  boot = FALSE,
  nboot = 10000,
  check.table = TRUE,
  ...
)
}
\arguments{
\item{x}{either an object of class "table", "matrix" or "numeric" representing a 2x2 contingency table,
or a "character" (a set of gene identifiers) or "list" or "tableList" object.
See the details section for more information.}

\item{...}{additional arguments for function \code{buildEnrichTable}.}

\item{dis}{Sorensen-Dice dissimilarity value. Only required to speed computations if this value
is known in advance.}

\item{se}{standard error estimate of the sample dissimilarity. Only required to speed computations
if this value is known in advance.}

\item{conf.level}{confidence level of the one-sided confidence interval, a numeric value between 0 and 1.}

\item{z.conf.level}{standard normal (or bootstrap, see arguments below) distribution quantile at the
\code{1 - conf.level} value. Only required to speed computations if this value is known in advance.
Then, the argument \code{conf.level} is ignored.}

\item{boot}{boolean. If TRUE, \code{z.conf.level} is computed by means of a bootstrap
approach instead of the asymptotic normal approach. Defaults to FALSE.}

\item{nboot}{numeric, number of initially planned bootstrap replicates. Ignored if
\code{boot == FALSE}. Defaults to 10000.}

\item{check.table}{Boolean. If TRUE (default), argument \code{x} is checked to adequately
represent a 2x2 contingency table. This checking is performed by means of function
\code{nice2x2Table}.}

\item{y}{an object of class "character" representing a vector of gene identifiers (e.g., ENTREZ).}
}
\value{
In the "table", "matrix", "numeric" and "character" interfaces, the value of the Upper
limit of the confidence interval for the Sorensen-Dice dissimilarity.
When \code{boot == TRUE}, this result also haves a an extra attribute: "eff.nboot" which
corresponds to the number of effective bootstrap replicats, see the details section.
In the "list" and "tableList" interfaces, the result is the symmetric matrix
of all pairwise upper limits.
}
\description{
Upper limit of a one-sided confidence interval (0, dUpp] for the Sorensen-Dice dissimilarity
}
\details{
This function computes the upper limit of a one-sided confidence interval for the Sorensen-Dice dissimilarity,
given a 2x2 arrangement of frequencies (either implemented as a "table", a "matrix"
or a "numeric" object):

\tabular{rr}{
\eqn{n_{11}} \tab \eqn{n_{10}} \cr
\eqn{n_{01}} \tab \eqn{n_{00}},
}

The subindex '11' corresponds to those
GO terms enriched in both lists, '01' to terms enriched in the second list but not in the first one,
'10' to terms enriched in the first list but not enriched in the second one and '00' corresponds
to those GO terms non enriched in both gene lists, i.e., to the double negatives, a value which
is ignored in the computations, except if \code{boot == TRUE}.

In the "numeric" interface, if \code{length(x) >= 4}, the values are interpreted
as
\eqn{(n_{11}, n_{01}, n_{10}, n_{00})}{%
(n_11, n_01, n_10, n_00)}, always in this order and discarding extra values if necessary.

Arguments \code{dis}, \code{se} and \code{z.conf.level} are not required. If known in advance (e.g., as
a consequence of previous computations with the same data), providing its value may speed the computations.

By default, \code{z.conf.level} corresponds to the 1 - conf.level quantile of a standard normal N(0,1)
distribution, as the studentized statistic (^d - d) / ^se) is asymptotically N(0,1). In
the studentized statistic, d stands for the "true" Sorensen-Dice dissimilarity, ^d to its sample estimate
and ^se for the estimate of its standard error.
In fact, the normal is its limiting distribution but, for finite samples, the true sampling
distribution may present departures from normality (mainly with some inflation in the
left tail).
The bootstrap method provides a better approximation to the true sampling distribution.
In the bootstrap approach, \code{nboot} new bootstrap contingency tables are generated from a
multinomial distribution with parameters
\code{size =} \eqn{n = n_{11} + n_{01} + n_{10} + n_{00}}{%
n11 + n01 + n10 + n00} and probabilities
\eqn{(n_{11} / n, n_{01} / n, n_{10}, n_{00} / n)}{\%}.
Sometimes, some of these generated tables may present so low
frequencies of enrichment that make them unable for Sorensen-Dice computations. As a consequence,
the number of effective bootstrap samples may be lower than the number of initially planned bootstrap
samples \code{nboot}.
Computing in advance the value of argument \code{z.conf.level} may be a way to cope with
these departures from normality, by means of a more adequate quantile function.
Alternatively, if \code{boot == TRUE}, a bootstrap quantile is internally computed.

If \code{x} is an object of class "character", then \code{x} (and \code{y}) must represent
two "character" vectors of valid gene identifiers (e.g., ENTREZ).
Then the confidence interval for the dissimilarity between lists \code{x} and \code{y} is computed,
after internally summarizing them as a 2x2 contingency table of joint enrichment.
This last operation is performed by function \code{\link{buildEnrichTable}} and "valid gene
identifiers (e.g., ENTREZ)" stands for the coherency of these gene identifiers with the arguments
\code{geneUniverse} and \code{orgPackg} of \code{buildEnrichTable}, passed by the ellipsis
argument \code{...} in \code{dUppSorensen}.

In the "list" interface, the argument must be a list of "character" vectors, each one
representing a gene list (character identifiers). Then, all pairwise upper limits of
the dissimilarity between these gene lists are computed.

In the "tableList" interface, the upper limits are computed over each one of these tables.
Given gene lists (i.e. "character" vectors of gene identifiers) l1, l2, ..., lk,
an object of class "tableList" (typically constructed by a call to function
\code{\link{buildEnrichTable}}) is a list of lists of
contingency tables t(i,j) generated from each pair of gene lists i and j, with the
following structure:

$l2

$l2$l1$t(2,1)

$l3

$l3$l1$t(3,1), $l3$l2$t(3,2)

...

$lk

$lk$l1$t(k,1), $lk$l2$t(k,2), ..., $lk$l(k-1)t(k,k-1)
}
\section{Methods (by class)}{
\itemize{
\item \code{duppSorensen(table)}: S3 method for class "table"

\item \code{duppSorensen(matrix)}: S3 method for class "matrix"

\item \code{duppSorensen(numeric)}: S3 method for class "numeric"

\item \code{duppSorensen(character)}: S3 method for class "character"

\item \code{duppSorensen(list)}: S3 method for class "list"

\item \code{duppSorensen(tableList)}: S3 method for class "tableList"

}}
\examples{
# Gene lists 'atlas' and 'sanger' in 'Cangenes' dataset. Table of joint enrichment
# of GO terms in ontology BP at level 3.
data(cont_atlas.sanger_BP4)
?cont_atlas.sanger_BP4
duppSorensen(cont_atlas.sanger_BP4)
dSorensen(cont_atlas.sanger_BP4) + qnorm(0.95) * seSorensen(cont_atlas.sanger_BP4)
# Using the bootstrap approximation instead of the normal approximation to
# the sampling distribution of (^d - d) / se(^d):
duppSorensen(cont_atlas.sanger_BP4, boot = TRUE)

# Contingency table as a numeric vector:
duppSorensen(c(56, 1, 30, 47))
duppSorensen(c(56, 1, 30))

# Upper confidence limit for the Sorensen-Dice dissimilarity, from scratch,
# directly from two gene lists:
# (These examples may be considerably time consuming due to many enrichment
# tests to build the contingency tables of mutual enrichment)
# data(allOncoGeneLists)
# ?allOncoGeneLists

# Obtaining ENTREZ identifiers for the gene universe of humans:
# library(org.Hs.eg.db)
# humanEntrezIDs <- keys(org.Hs.eg.db, keytype = "ENTREZID")

# Computing the Upper confidence limit:
# duppSorensen(allOncoGeneLists$atlas, allOncoGeneLists$sanger,
#              onto = "CC", GOLevel = 5,
#              geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
# Even more time consuming (all pairwise values):
# duppSorensen(allOncoGeneLists,
#              onto = "CC", GOLevel = 5,
#              geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")
}
\seealso{
\code{\link{buildEnrichTable}} for constructing contingency tables of mutual
enrichment,
\code{\link{nice2x2Table}} for checking contingency tables validity,
\code{\link{dSorensen}} for computing the Sorensen-Dice dissimilarity,
\code{\link{seSorensen}} for computing the standard error of the dissimilarity,
\code{\link{equivTestSorensen}} for an equivalence test.
}
