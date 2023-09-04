# goSorensen
## version 1.1.0
- This is the first version submitted to review. There are currently no updates

## version 1.2.0
- Added some results (allTabsBP.4 object) about computed crosstabs for all possible list pairs in allOncoGeneLists for BP ontology, level 4

## version 1.3.0
We add the following functions

- __allBuildEnrichTable__: Given k lists of genes, it generates the k(k–1)/2 contingency
tables of joint enrichment for all possible pairs of lists.
- __allEquivTestSorensen__: Accepts the objects created with allBuildEnrichTable as
its first argument and quickly obtains the k(k–1)/2 equivalence tests.
- __sorenThreshold__: Implements an algorithm allowing the computation of the “equivalence threshold” dissimilarities matrix for all the k(k–1)/2 tests. According to what is indicated in the arguments, the results of this function are stored in a matrix for a specific ontology and level or a list with a matrix for more than one ontology and/or level.
- __hclustThreshold__: Generates an object of class "hclust". For a specific ontology
and level, plots a dendrogram where all the k(k–1)/2 comparisons are joined at
the height of their respective "equivalence threshold" dissimilarity.
- __allHclustThreshold__: Performs the same calculations as hclustThreshold but for the specified GO ontologies and levels
(all three ontologies PB, MF, and CC and levels from 2 to 10 by default)

In addition: 

- We improve the vignette "An introduction to the goSorensen package" by implementing the new functions mentioned above and updating the results of the examples with the results of the latest version of Bioconductor.
