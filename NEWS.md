# goSorensen
## version 1.1.0
- This is the first version submitted to review. There are currently no updates.

## version 1.2.0
- Added some results (allTabsBP.4 object) about computed crosstabs for all possible list pairs in allOncoGeneLists for BP ontology, level 4

## version 1.3.0
We add the following functions:

- __allBuildEnrichTable__: This function produces k(k–1)/2 contingency tables by comparing all possible pairs of feature lists using the provided GO ontologies and GO levels.
- __allEquivTestSorensen__: This function computes the k(k–1)/2 equivalence tests by comparing all possible pairs of feature lists using the provided GO ontologies and GO levels.
- __sorenThreshold__: This function applies an algorithm to compute the irrelevance-threshold matrix of dissimilarities for a specific ontology and GO level.
- __allSorenThreshold__: This function iterates __sorenThreshold__ along the specified GO ontologies and GO levels
- __hclustThreshold__: This function plots a dendrogram for a specific ontology and GO level based on the irrelevance-threshold matrix of dissimilarities..
- __allHclustThreshold__: This function iterates __hclustThreshold__ along the specified GO ontologies and GO levels

In addition: 

- We enhance the vignette "An introduction to the goSorensen package" by incorporating the newly introduced functions and updating the example results with the latest version of Bioconductor. 

## version 1.8.0
We add the following function:

- __enrichedIn__: This function generates a binary matrix in which the rows represent GO terms, and the columns represent lists. The matrix uses the values TRUE or FALSE to indicate whether each GO term is enriched or not for each list. 

In addition: 

- We have updated the __buildEnrichedTable__ and __allBuildEnrichedTable__ functions. They now uses the outcomes of __enrichedIn__ to build the contingency table. During the conducted tests, the speed of the procedure for generating contingency tables increased by a factor of six compared to the previous version. 

  Furthermore, we have included the _showEnrichedIn_ argument in these functions. This is a boolean argument. If the argument is TRUE, in addition to the enrichment contingency table, the function saves a matrix in the global environment, which contains the cross table of the enriched and non-enriched GO terms vs the names of the gene lists generated with the __enrichedIn__ function. This matrix is stored under the name of "enrichedIn_" followed by the name of the ontology and the level being analysed. An object will be produced for every scenario if there are several levels and ontologies. 

- We add the vignette "__irrelevance-threshold_Matrix_Dissimilarities__." This vignette illustrates calculating, visualizing, and interpreting the irrelevance-threshold matrix of dissimilarities _D_. This matrix provides dissimilarities between pairs of compared lists. These dissimilarities are not only a descriptive measure but also based on the irrelevance threshold determining whether two lists are equivalent. So, this dissimilarity measure between the two lists is directly associated with their declaration of equivalence.

  The matrix _D_ can be represented in interpretable statistic graphs such as dendrograms or biplots, which help to visualize the formation of groups containing equivalent lists. Furthermore, interpreting the biplot dimensions gives us the biological functions responsible for the equivalence between lists.

## version 1.10.0
In this updated version, unless the user specifies otherwise in the _onlyEnriched_ argument of the _enrichedIn_ function, the enrichment matrix of GO terms includes only those terms that are enriched in at least one of the lists being compared, excluding all terms that are not enriched in any of the lists. This optimization significantly reduces computation times for enrichment analysis compared to previous versions, such as when generating enrichment contingency tables using _buildEnrichTable_.

Additionally, with the latest versions of Bioconductor and its updated packages, we have revised the naming conventions and outputs for results obtained from the primary __goSorensen__ functions. These changes are summarized in two new bullet points, offering clearer guidance and improved usability for package users compared to earlier versions.