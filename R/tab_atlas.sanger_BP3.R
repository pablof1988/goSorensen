#' Cross-tabulation of enriched GO items at level 3 of ontology BP in two gene lists
#'
#' From the "Cancer gene list" of Bushman Lab, a collection of gene lists related with cancer,
#' for gene lists "Atlas" and "Sanger", this dataset is the cross-tabulation of all GO items
#' of ontology BP at level 3 which are:
#' non-enriched in both lists, enriched in atlas but not in sanger, non-enriched in atlas
#' but enriched in sanger and enriched in both lists.
#'
#' @format An object of class "table" representing a 2x2 contingency table.
#' @source \url{http://www.bushmanlab.org/links/genelists}
"tab_atlas.sanger_BP3"
