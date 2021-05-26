adjSignifPvals <- function(x, alpha = 0.05, adj.meth = "holm") {
  rawPvals <- getPvalue(x)
  ontoNams <- names(rawPvals)
  adjPvals <- lapply(ontoNams, function(onto){
    levNams <- names(rawPvals[[onto]])
    adjLevPvals <- lapply(levNams, function(lev){
      pvals <- rawPvals[[onto]][[lev]]
      pvals <- p.adjust(pvals[!is.na(pvals)], method = adj.meth)
      pvals <- pvals[pvals <= alpha]
      return(pvals)
    })
    names(adjLevPvals) <- levNams
    return(adjLevPvals)
  })
  names(adjPvals) <- ontoNams

  result <- lapply(ontoNams, function(onto){
    levNams <- names(x[[onto]])
    levResult <- lapply(levNams, function(lev){
      signifPairNams <- names(adjPvals[[onto]][[lev]])
      signifPairs <- vector(mode = "list", length = length(signifPairNams))
      names(signifPairs) <- signifPairNams
      for (iLst1 in names(x[[onto]][[lev]])){
        for (iLst2 in names(x[[onto]][[lev]][[iLst1]])) {
          thisPair <- paste(iLst1, iLst2, sep = ".")
          if (thisPair %in% signifPairNams) {
            signifPairs[[thisPair]] <- list(p.value = adjPvals[[onto]][[lev]][thisPair],
                                            enrichTab = getTable(x[[onto]][[lev]][[iLst1]][[iLst2]]))
          }
        }
      }
      return(signifPairs)
    })
    names(levResult) <- levNams
    return(levResult)
  })
  names(result) <- ontoNams
  return(result)
}
