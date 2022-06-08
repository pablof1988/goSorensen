# goSorensen
This R package implements inferential methods to compare gene lists (in this first release, to prove equivalence) in terms of their biological meaning as expressed in the GO. The compared gene lists are characterized by cross-tabulation frequency tables of enriched GO items. Dissimilarity between gene lists is evaluated using the Sorensen-Dice index.
The fundamental guiding principle is that two gene lists are taken as similar if they share a great proportion of common enriched GO items.

This inferential method, is developed and explained in the paper *An equivalence test between features lists, based on the Sorensen-Dice index and the joint frequencies of GO term enrichment*, available online in <https://rdcu.be/cOISz>

## Installation instructions
goSorensen package has to be installed with a working R version (>=4.2.0). Installation could take a few minutes on a regular desktop or laptop. Package can be installed from Bioconductor or `devtools` package, then it needs to be loaded using `library(goSorensen)`

To install from Bioconductor (recommended):
```{r}
## Only if BiocManager is not previosly installed:
install.packages("BiocManager")

## otherwise, directly:
BiocManager::install("goSorensen")
```

To install from Github

```{r}
devtools::install_github("pablof1988/goSorensen", build_vignettes = TRUE)
```

## Basic ideas and concepts
For the moment, the goSorensen package provides the following functions:

- **buildEnrichTable:** Build an enrichment contingency table from two gene lists. 
- **nice2x2Table:** Check for validity an enrichment contingency table.
- **dSorensen:** Compute the Sorensen-Dice dissimilarity.
- **seSorensen:** Standard error estimate of the sample Sorensen-Dice dissimilarity.
- **duppSorensen:** Upper limit of a one-sided confidence interval (0,dUpp] for the population dissimilarity.
- **equivTestSorensen:** Equivalence test between two gene lists, based on the Sorensen-Dice dissimilarity.
- **allEquivTestSorensen** Iterate equivTestSorensen along GO ontologies and GO levels.
- **getDissimilarity, getPvalue, getSE, getTable, getUpper, getNboot:** Accessor functions to some fields of an equivalence test result
- **upgrade:** Updating the result of an equivalence test, e.g., changing the equivalence limit.


All these functions are generic, with methods for classes representing diverse kinds of data or statistical results. dSorensen, seSorensen, duppSorensen and equivTestSorensen have methods to manage classes directly representing enrichment cross-tabulations or, alternatively, character vectors representing gene lists. In this second case (gene lists), the enrichment tables are built in a previous step. If the first parameter of these functions is a list of character vectors (i.e., a list of gene lists) all paired enrichment cross-tabulations are performded and the result is a symmetric matrix of all paired dissimilarities, or standard errors, or confidence interval upper limits -or a list structure emulating a symmetric matrix as a result of equivTestSorensen. The accessor and upgradding functions have methods for all classes of equivTestSorensen and allEquivTestSorensen results.

## Contribution Guidelines
Contributions are welcome, if you wish to contribute or give ideas to improve the package, please you can contact with maintainer (Pablo Flores) to the addres `p_flores@espoch.edu.ec`, and we can discuss your suggestion.

## References
<div id="refs" class="references">
<div id="goSorensen">

Flores, P., Salicru, M., Sanchez-Pla, A., & Ocana, J. (2022). An equivalence test between features lists, based on the Sorensen-Dice index and the joint frequencies of GO term enrichment. BMC bioinformatics, 23(1), 1-21.

</div>
</div>
