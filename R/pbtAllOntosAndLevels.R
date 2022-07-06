#' The Sorensen-Dice test performed on some gene lists possibly related to kidney rejection after transplantation
#' based on non-updated information at \url{https://rdrr.io/cran/tcgsaseq/man/PBT_gmt.html}, take them just as an
#' illustrative example.
#' Tests performed for all GO ontologies and for GO levels 3 to 10.
#'
#' For each ontology and GO level, the result contains the result of all pairwise tests of equivalence between
#' these gene lists.
#'
#' @format An object of class "AllEquivSDhtest" inheriting from class "list". Each one of its elements, named
#' BP, CC and MF respectively, corresponds to a GO ontology. It is itself a list of length 8 whose elements
#' are named as "Level 3" to "Level 10". For each combination of ontology and level, there is an object of
#' class "equivSDhtestList" codifying the result of all pairwise tests between these kidney rejection gene
#' lists.
#'
#' @examples
#' # This code may help to understand the structure of these data:
#' data(pbtAllOntosAndLevels)
#' ?pbtAllOntosAndLevels
#' names(pbtAllOntosAndLevels)
#' names(pbtAllOntosAndLevels$BP)
#' names(pbtAllOntosAndLevels$BP$`level 4`)
#' class(pbtAllOntosAndLevels$BP$`level 4`)
#' pbtAllOntosAndLevels$BP$`level 4`
#' names(pbtAllOntosAndLevels$BP$`level 4`$KT1)
#' class(pbtAllOntosAndLevels$BP$`level 4`$KT1)
#' class(pbtAllOntosAndLevels$BP$`level 4`$KT1$IRITD5)
#' pbtAllOntosAndLevels$BP$`level 4`$KT1$IRITD5
#' @usage data(pbtAllOntosAndLevels)
"pbtAllOntosAndLevels"
