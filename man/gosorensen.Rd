% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gosorensen.R
\docType{package}
\name{gosorensen}
\alias{goSorensen}
\alias{goSorensen-package}
\alias{gosorensen}
\title{gosorensen: A package for making inference on gene lists based on the Sorensen-Dice dissimilarity}
\description{
Given two lists of genes, and a set of Gene Ontology (GO) items (e.g., all GO items in a given
level of a given GO ontology) one may explore some aspects of their biological meaning by
constructing a 2x2 contingency table, the cross-tabulation of: number of these GO items
non-enriched in both gene lists (n00), items enriched in the first list but not in the second one
(n10), items non-enriched in the first list but enriched in the second (n10) and items enriched
in both lists (n11). Then, one may express the degree of similarity or dissimilarity between the
two lists by means of an appropriate index computed on these frequency tables of concordance or
non-concordance in GO items enrichment.
In our opinion, an appropriate index is the Sorensen-Dice index which ignores the double negatives
n00: if the total number of candidate GO items under consideration grows (e.g., all items in a deep
level of an ontology) likely n00 will also grow artificially. On the other hand, intuitively the
degree of similarity between both lists must be directly related to the degree of concordance in
the enrichment, n11.
}
\details{
gosorensen package provides the following functions:
\describe{
  \item{enrichedIn}{Build a cross-tabulation of enriched and non-enriched GO terms vs. gene lists}
  \item{buildEnrichTable}{Build an enrichment contingency table from two or more gene lists}
  \item{allBuildEnrichTable}{Iterate 'buildEnrichTable' along the specified GO ontologies and GO levels}
  \item{nice2x2Table}{Check for validity an enrichment contingency table}
  \item{dSorensen}{Compute the Sorensen-Dice dissimilarity}
  \item{seSorensen}{Standard error estimate of the sample Sorensen-Dice dissimilarity}
  \item{duppSorensen}{Upper limit of a one-sided confidence interval (0,dUpp] for the population dissimilarity}
  \item{equivTestSorensen}{Equivalence test between two gene lists, based on the Sorensen-Dice dissimilarity}
  \item{allEquivTestSorensen}{Iterate equivTestSorensen along GO ontologies and GO levels}
  \item{getDissimilarity, getPvalue, getSE, getTable, getUpper, getNboot, getEffNboot}{Accessor functions to some fields of an equivalence test result}
  \item{upgrade}{Updating the result of an equivalence test, e.g., changing the equivalence limit}
  \item{sorenThreshold}{For a given level (2, 3, ...) in a GO ontology (BP, MF or CC), compute the equivalence threshold dissimilarity matrix.}
  \item{allSorenThreshold}{Iterate 'sorenThreshold' along the specified GO ontologies and GO levels.}
  \item{hclustThreshold}{From a Sorensen-Dice threshold dissimilarity matrix, generate an object of class "hclust"}
  \item{allHclustThreshold}{Iterate 'hclustThreshold' along the specified GO ontologies and GO levels}
  \item{pruneClusts}{Remove all NULL or unrepresentable as a dendrogram "equivClustSorensen" elements in an object of class "equivClustSorensenList"}
}
All these functions are generic, adequate for different (S3) classes representing the before cited
GO term enrichment cross-tabulations.
}
\author{
\strong{Maintainer}: Pablo Flores \email{p_flores@espoch.edu.ec} (\href{https://orcid.org/0000-0002-7156-8547}{ORCID})

Authors:
\itemize{
  \item Jordi Ocana (\href{https://orcid.org/0000-0002-4736-699}{ORCID}) [contributor]
}

Other contributors:
\itemize{
  \item Alexandre Sanchez-Pla (\href{https://orcid.org/0000-0002-8673-7737}{ORCID}) [contributor]
  \item Miquel Salicru (\href{https://orcid.org/0000-0001-9644-5626}{ORCID}) [contributor]
}

}
