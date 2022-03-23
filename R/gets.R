#' Access to specific fields of an equivalence test result
#'
#' Accessor functions to specific fields of objects of classes "equivSDhtest", "equivSDhtestList"
#' or "allEquivSDtest", i.e., to the result of functions 'equivTestSorensen' and 'allEquivTestSorensen'.
#'
#' @param x an object of class "equivSDhtest" or "equivSDhtestList" or "allEquivSDtest".
#' @param onto character, a vector with one or more of "BP", "CC" or "MF", ontologies to access.
#' @param GOLevel numeric, a vector with one or more GO levels to access.
#' @param simplify logical, if TRUE the result is simplified, e.g., returning a vector instead
#' of a matrix.
#' @param listNames character(2), the names of a pair of gene lists.
#' @examples
#' # Result of the equivalence test between gene lists 'waldman' and 'atlas', in dataset
#' # 'allOncoGeneLists', at level 4 of the BP ontology:
#' waldman_atlas.BP.4
#' class(waldman_atlas.BP.4)
#' # This may correspond to the result of code like:
#' # waldman_atlas.BP.4 <- equivTestSorensen(
#' #   allOncoGeneLists[["waldman"]], allOncoGeneLists[["atlas"]],
#' #   geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#' #   onto = "BP", GOLevel = 4, listNames = c("waldman", "atlas"))
#' getPvalue(waldman_atlas.BP.4)
#' getDissimilarity(waldman_atlas.BP.4)
#' getUpper(waldman_atlas.BP.4)
#' getSE(waldman_atlas.BP.4)
#' getTable(waldman_atlas.BP.4)
#'
#' # All pairwise equivalence tests at level 4 of the BP ontology
#' BP.4
#' class(BP.4)
#' # This may correspond to a call like:
#' # BP.4 <- equivTestSorensen(allOncoGeneLists,
#' #                           geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db",
#' #                           onto = "BP", GOLevel = 4)
#' getPvalue(BP.4)
#' getPvalue(BP.4, simplify = FALSE)
#' getDissimilarity(BP.4)
#' getDissimilarity(BP.4, simplify = FALSE)
#' getUpper(BP.4)
#' getUpper(BP.4, simplify = FALSE)
#' getSE(BP.4)
#' getSE(BP.4, simplify = FALSE)
#' getTable(BP.4)
#'
#' # Equivalence test iteradted over all GO ontologies and levels 3 to 10.
#' cancerEquivSorensen
#' class(cancerEquivSorensen)
#' # This may correspond to code like:
#' # (By default, the tests are iterated over all GO ontologies and for levels 3 to 10)
#' # cancerEquivSorensen <- allEquivTestSorensen(allOncoGeneLists,
#' #                                             geneUniverse = humanEntrezIDs,
#' #                                             orgPackg = "org.Hs.eg.db")

#' # 2x2 contingecy tables of joint enrichment:
#' getTable(cancerEquivSorensen)
#' getTable(cancerEquivSorensen, GOLevel = "level 6")
#' getTable(cancerEquivSorensen, GOLevel = "level 6", listNames = c("waldman", "sanger"))
#' getTable(cancerEquivSorensen, GOLevel = "level 6", onto = "BP")
#' getTable(cancerEquivSorensen, GOLevel = "level 6", onto = "BP",
#'          listNames = c("waldman", "sanger"))
#'
#' # p-values:
#' getPvalue(cancerEquivSorensen)
#' getPvalue(cancerEquivSorensen, simplify = FALSE)
#' getPvalue(cancerEquivSorensen, GOLevel = "level 6")
#' getPvalue(cancerEquivSorensen, GOLevel = "level 6", simplify = FALSE)
#' getPvalue(cancerEquivSorensen, GOLevel = "level 6", listNames = c("waldman", "sanger"))
#' getPvalue(cancerEquivSorensen, GOLevel = "level 6", onto = "BP")
#' getPvalue(cancerEquivSorensen, GOLevel = "level 6", onto = "BP", simplify = FALSE)
#' getPvalue(cancerEquivSorensen, GOLevel = "level 6", onto = "BP",
#'           listNames = c("waldman", "sanger"))
#' getPvalue(cancerEquivSorensen$BP$`level 4`)
#'
#' # Estimated Sorensen-Dice dissimilarity:
#' getDissimilarity(cancerEquivSorensen)
#' getDissimilarity(cancerEquivSorensen, simplify = FALSE)
#' getDissimilarity(cancerEquivSorensen, GOLevel = "level 6")
#' getDissimilarity(cancerEquivSorensen, GOLevel = "level 6", simplify = FALSE)
#' getDissimilarity(cancerEquivSorensen, GOLevel = "level 6", listNames = c("waldman", "sanger"))
#' getDissimilarity(cancerEquivSorensen, GOLevel = "level 6", onto = "BP")
#' getDissimilarity(cancerEquivSorensen, GOLevel = "level 6", onto = "BP", simplify = FALSE)
#' getDissimilarity(cancerEquivSorensen, GOLevel = "level 6", onto = "BP",
#'                  listNames = c("waldman", "sanger"))
#' getDissimilarity(cancerEquivSorensen$BP$`level 4`)
#'
#' # Upper confidence limits for the Sorensen-Dice dissimilarity:
#' getUpper(cancerEquivSorensen)
#' getUpper(cancerEquivSorensen, simplify = FALSE)
#' getUpper(cancerEquivSorensen, GOLevel = "level 6")
#' getUpper(cancerEquivSorensen, GOLevel = "level 6", simplify = FALSE)
#' getUpper(cancerEquivSorensen, GOLevel = "level 6", listNames = c("waldman", "sanger"))
#' getUpper(cancerEquivSorensen, GOLevel = "level 6", onto = "BP")
#' getUpper(cancerEquivSorensen, GOLevel = "level 6", onto = "BP", simplify = FALSE)
#' getUpper(cancerEquivSorensen, GOLevel = "level 6", onto = "BP",
#'          listNames = c("waldman", "sanger"))
#' getUpper(cancerEquivSorensen$BP$`level 4`)
#'
#' # Standard error of the Sorensen-Dice dissimilarity estimate:
#' getSE(cancerEquivSorensen)
#' getSE(cancerEquivSorensen, simplify = FALSE)
#' getSE(cancerEquivSorensen, GOLevel = "level 6")
#' getSE(cancerEquivSorensen, GOLevel = "level 6", simplify = FALSE)
#' getSE(cancerEquivSorensen, GOLevel = "level 6", listNames = c("waldman", "sanger"))
#' getSE(cancerEquivSorensen, GOLevel = "level 6", onto = "BP")
#' getSE(cancerEquivSorensen, GOLevel = "level 6", onto = "BP", simplify = FALSE)
#' getSE(cancerEquivSorensen, GOLevel = "level 6", onto = "BP", listNames = c("waldman", "sanger"))
#' getSE(cancerEquivSorensen$BP$`level 4`)
#'
#'
#' # Acces to Sorensen-Dice asymptotic test enrichment 2x2 contingency table
#'
#' @export
getTable <- function(x, ...) {
  UseMethod("getTable")
}

#' @describeIn getTable S3 method for class "equivSDhtest"
#' @export
getTable.equivSDhtest <- function(x) {
  return(x$enrichTab)
}

#' @describeIn getTable S3 method for class "equivSDhtestList"
#' @export
getTable.equivSDhtestList <- function(x) {
  result <- lapply(x, function(xi){
    resaux <- lapply(xi, getTable.equivSDhtest)
    names(resaux) <- names(xi)
    return(resaux)
  })
  names(result) <- names(x)
  return(result)
  # return(lapply(x, getTable.equivSDhtest))
}

#' @describeIn getTable S3 method for class "AllEquivSDhtest"
#' @export
getTable.AllEquivSDhtest <- function(x,
                                     onto, GOLevel, listNames) {
  if (missing(onto)) {
    onto <- names(x)
  }
  if (missing(GOLevel)) {
    GOLevel <- names(x[[1]])
  }
  allLists <- missing(listNames)
  result <- lapply(onto, function(ionto) {
    resLev <- lapply(GOLevel, function(ilev) {
      if (allLists) {
        namsList1 <- names(x[[ionto]][[ilev]])
        resList1 <- lapply(namsList1, function(ilist1) {
          namsList2 <- names(x[[ionto]][[ilev]][[ilist1]])
          resList2 <- lapply(namsList2, function(ilist2) {
            return(x[[ionto]][[ilev]][[ilist1]][[ilist2]]$enrichTab)
          })
          names(resList2) <- namsList2
          return(resList2)
        })
        names(resList1) <- namsList1
        # if (simplify) {
        #   resList1 <- unlist(resList1, recursive = FALSE)
        # }
        return(resList1)
      }else {
        return(x[[ionto]][[ilev]][[listNames[1]]][[listNames[2]]]$enrichTab)
      }
    })
    names(resLev) <- GOLevel
    # if (simplify && allLists) {
    #   resLev <- unlist(resLev, recursive = FALSE)
    # }
    return(resLev)
  })
  names(result) <- onto
  # if (simplify) {
  #   result <- unlist(result, recursive = FALSE)
  # }
  return(result)
}

#' Acces to Sorensen-Dice equivalence test p-value
#'
#' @export
getPvalue <- function(x, ...) {
  UseMethod("getPvalue")
}

#' @describeIn getPvalue S3 method for class "equivSDhtest"
#' @export
getPvalue.equivSDhtest <- function(x) {
  return(x$p.value[1])
}

#' @describeIn getPvalue S3 method for class "equivSDhtestList"
#' @export
getPvalue.equivSDhtestList <- function(x, simplify = TRUE) {
  result <- lapply(x, function(xi){
    resaux <- lapply(xi, getPvalue.equivSDhtest)
    names(resaux) <- names(xi)
    return(resaux)
  })
  namsList1 <- names(result) <- names(x)
  result <- unlist(result)
  if (!simplify) {
    namsMat <- c(names(x[[1]])[1], namsList1)
    nrows <- length(namsList1) + 1
    resMat <- matrix(0.0, nrow = nrows, ncol = nrows)
    resMat[upper.tri(resMat)] <- result
    resMat[lower.tri(resMat)] <- t(resMat)[lower.tri(resMat)]
    rownames(resMat) <- namsMat
    colnames(resMat) <- namsMat
    return(resMat)
  }else{
    return(result)
  }
}

#' @describeIn getPvalue S3 method for class "AllEquivSDhtest"
#' @export
getPvalue.AllEquivSDhtest <- function(x,
                                      onto, GOLevel, listNames,
                                      simplify = TRUE) {
  if (missing(onto)) {
    onto <- names(x)
  }
  if (missing(GOLevel)) {
    GOLevel <- names(x[[1]])
  }
  allLists <- missing(listNames)
  result <- lapply(onto, function(ionto) {
    resLev <- lapply(GOLevel, function(ilev) {
      if (allLists) {
        namsList1 <- names(x[[ionto]][[ilev]])
        namsMat <- c(names(x[[ionto]][[ilev]][[1]])[1], namsList1)
        resList1 <- sapply(namsList1, function(ilist1) {
          namsList2 <- names(x[[ionto]][[ilev]][[ilist1]])
          resList2 <- sapply(namsList2, function(ilist2) {
            return(x[[ionto]][[ilev]][[ilist1]][[ilist2]]$p.value)
          })
          names(resList2) <- namsList2
          return(resList2)
        })
        names(resList1) <- namsList1
        resList1 <- unlist(resList1, recursive = FALSE)
        if (!simplify) {
          nrows <- length(namsList1) + 1
          matList1 <- matrix(0.0, nrow = nrows, ncol = nrows)
          matList1[upper.tri(matList1)] <- resList1
          matList1[lower.tri(matList1)] <- t(matList1)[lower.tri(matList1)]
          rownames(matList1) <- namsMat
          colnames(matList1) <- namsMat
          resList1 <- matList1
        }
        return(resList1)
      }else {
        return(x[[ionto]][[ilev]][[listNames[1]]][[listNames[2]]]$p.value)
      }
    })
    names(resLev) <- GOLevel
    # if (simplify) {
    #   resLev <- unlist(resLev, recursive = FALSE)
    # }
    return(resLev)
  })
  names(result) <- onto
  # if (simplify) {
  #   result <- unlist(result, recursive = FALSE)
  # }
  return(result)
}

#' Acces to the Sorensen-Dice dissimilarity
#'
#' @export
getDissimilarity <- function(x, ...) {
  UseMethod("getDissimilarity")
}

#' @describeIn getDissimilarity S3 method for class "equivSDhtest"
#' @export
getDissimilarity.equivSDhtest <- function(x) {
  return(x$estimate)
}

#' @describeIn getDissimilarity S3 method for class "equivSDhtestList"
#' @export
getDissimilarity.equivSDhtestList <- function(x, simplify = TRUE) {
  result <- lapply(x, function(xi){
    resaux <- lapply(xi, getDissimilarity.equivSDhtest)
    names(resaux) <- names(xi)
    return(resaux)
  })
  namsList1 <- names(result) <- names(x)
  result <- unlist(result)
  if (!simplify) {
    namsMat <- c(names(x[[1]])[1], namsList1)
    nrows <- length(namsList1) + 1
    resMat <- matrix(0.0, nrow = nrows, ncol = nrows)
    resMat[upper.tri(resMat)] <- result
    resMat[lower.tri(resMat)] <- t(resMat)[lower.tri(resMat)]
    rownames(resMat) <- namsMat
    colnames(resMat) <- namsMat
    return(resMat)
  }else{
    return(result)
  }
}
#'
#' @describeIn getDissimilarity S3 method for class "AllEquivSDhtest"
#' @export
getDissimilarity.AllEquivSDhtest <- function(x, onto, GOLevel, listNames,
                                             simplify = TRUE) {
  if (missing(onto)) {
    onto <- names(x)
  }
  if (missing(GOLevel)) {
    GOLevel <- names(x[[1]])
  }
  allLists <- missing(listNames)
  result <- lapply(onto, function(ionto) {
    resLev <- lapply(GOLevel, function(ilev) {
      if (allLists) {
        namsList1 <- names(x[[ionto]][[ilev]])
        namsMat <- c(names(x[[ionto]][[ilev]][[1]])[1], namsList1)
        resList1 <- sapply(namsList1, function(ilist1) {
          namsList2 <- names(x[[ionto]][[ilev]][[ilist1]])
          resList2 <- sapply(namsList2, function(ilist2) {
            return(x[[ionto]][[ilev]][[ilist1]][[ilist2]]$estimate)
          })
          names(resList2) <- namsList2
          return(resList2)
        })
        names(resList1) <- namsList1
        resList1 <- unlist(resList1, recursive = FALSE)
        if (!simplify) {
          nrows <- length(namsList1) + 1
          matList1 <- matrix(0.0, nrow = nrows, ncol = nrows)
          matList1[upper.tri(matList1)] <- resList1
          matList1[lower.tri(matList1)] <- t(matList1)[lower.tri(matList1)]
          rownames(matList1) <- namsMat
          colnames(matList1) <- namsMat
          resList1 <- matList1
        }
        return(resList1)
      }else {
        return(x[[ionto]][[ilev]][[listNames[1]]][[listNames[2]]]$estimate)
      }
    })
    names(resLev) <- GOLevel
    # if (simplify) {
    #   resLev <- unlist(resLev, recursive = FALSE)
    # }
    return(resLev)
  })
  names(result) <- onto
  # if (simplify) {
  #   result <- unlist(result, recursive = FALSE)
  # }
  return(result)
}

#' Acces to Sorensen-Dice confidence interval upper limit
#'
#' @export
getUpper <- function(x, ...) {
  UseMethod("getUpper")
}

#' @describeIn getUpper S3 method for class "equivSDhtest"
#' @export
getUpper.equivSDhtest <- function(x) {
  return(x$conf.int[2])
}

#' @describeIn getUpper S3 method for class "equivSDhtestList"
#' @export
getUpper.equivSDhtestList <- function(x, simplify = TRUE) {
  result <- lapply(x, function(xi){
    resaux <- lapply(xi, getUpper.equivSDhtest)
    names(resaux) <- names(xi)
    return(resaux)
  })
  namsList1 <- names(result) <- names(x)
  result <- unlist(result)
  if (!simplify) {
    namsMat <- c(names(x[[1]])[1], namsList1)
    nrows <- length(namsList1) + 1
    resMat <- matrix(0.0, nrow = nrows, ncol = nrows)
    resMat[upper.tri(resMat)] <- result
    resMat[lower.tri(resMat)] <- t(resMat)[lower.tri(resMat)]
    rownames(resMat) <- namsMat
    colnames(resMat) <- namsMat
    return(resMat)
  }else{
    return(result)
  }
}
#'
#' @describeIn getUpper S3 method for class "AllEquivSDhtest"
#' @export
getUpper.AllEquivSDhtest <- function(x, onto, GOLevel, listNames,
                                     simplify = TRUE) {
  if (missing(onto)) {
    onto <- names(x)
  }
  if (missing(GOLevel)) {
    GOLevel <- names(x[[1]])
  }
  allLists <- missing(listNames)
  result <- lapply(onto, function(ionto) {
    resLev <- lapply(GOLevel, function(ilev) {
      if (allLists) {
        namsList1 <- names(x[[ionto]][[ilev]])
        namsMat <- c(names(x[[ionto]][[ilev]][[1]])[1], namsList1)
        resList1 <- sapply(namsList1, function(ilist1) {
          namsList2 <- names(x[[ionto]][[ilev]][[ilist1]])
          resList2 <- sapply(namsList2, function(ilist2) {
            return(x[[ionto]][[ilev]][[ilist1]][[ilist2]]$conf.int[2])
          })
          names(resList2) <- namsList2
          return(resList2)
        })
        names(resList1) <- namsList1
        resList1 <- unlist(resList1, recursive = FALSE)
        if (!simplify) {
          nrows <- length(namsList1) + 1
          matList1 <- matrix(0.0, nrow = nrows, ncol = nrows)
          matList1[upper.tri(matList1)] <- resList1
          matList1[lower.tri(matList1)] <- t(matList1)[lower.tri(matList1)]
          rownames(matList1) <- namsMat
          colnames(matList1) <- namsMat
          resList1 <- matList1
        }
        return(resList1)
      }else {
        return(x[[ionto]][[ilev]][[listNames[1]]][[listNames[2]]]$conf.int[2])
      }
    })
    names(resLev) <- GOLevel
    # if (simplify) {
    #   resLev <- unlist(resLev, recursive = FALSE)
    # }
    return(resLev)
  })
  names(result) <- onto
  # if (simplify) {
  #   result <- unlist(result, recursive = FALSE)
  # }
  return(result)
}


#' Acces to Sorensen-Dice asymptotic test standard error of the dissimilarity estimate
#'
#' @export
getSE <- function(x, ...) {
  UseMethod("getSE")
}

#' @describeIn getSE S3 method for class "equivSDhtest"
#' @export
getSE.equivSDhtest <- function(x) {
  return(attr(x$estimate, "se"))
}

#' @describeIn getSE S3 method for class "equivSDhtestList"
#' @export
getSE.equivSDhtestList <- function(x, simplify = TRUE) {
  result <- lapply(x, function(xi){
    resaux <- lapply(xi, getSE.equivSDhtest)
    names(resaux) <- names(xi)
    return(resaux)
  })
  namsList1 <- names(result) <- names(x)
  result <- unlist(result)
  if (!simplify) {
    namsMat <- c(names(x[[1]])[1], namsList1)
    nrows <- length(namsList1) + 1
    resMat <- matrix(0.0, nrow = nrows, ncol = nrows)
    resMat[upper.tri(resMat)] <- result
    resMat[lower.tri(resMat)] <- t(resMat)[lower.tri(resMat)]
    rownames(resMat) <- namsMat
    colnames(resMat) <- namsMat
    return(resMat)
  }else{
    return(result)
  }
}
#'
#' @describeIn getSE S3 method for class "AllEquivSDhtest"
#' @export
getSE.AllEquivSDhtest <- function(x, onto, GOLevel, listNames,
                                     simplify = TRUE) {
  if (missing(onto)) {
    onto <- names(x)
  }
  if (missing(GOLevel)) {
    GOLevel <- names(x[[1]])
  }
  allLists <- missing(listNames)
  result <- lapply(onto, function(ionto) {
    resLev <- lapply(GOLevel, function(ilev) {
      if (allLists) {
        namsList1 <- names(x[[ionto]][[ilev]])
        namsMat <- c(names(x[[ionto]][[ilev]][[1]])[1], namsList1)
        resList1 <- sapply(namsList1, function(ilist1) {
          namsList2 <- names(x[[ionto]][[ilev]][[ilist1]])
          resList2 <- sapply(namsList2, function(ilist2) {
            return(attr(x[[ionto]][[ilev]][[ilist1]][[ilist2]]$estimate, "se"))
          })
          names(resList2) <- namsList2
          return(resList2)
        })
        names(resList1) <- namsList1
        resList1 <- unlist(resList1, recursive = FALSE)
        if (!simplify) {
          nrows <- length(namsList1) + 1
          matList1 <- matrix(0.0, nrow = nrows, ncol = nrows)
          matList1[upper.tri(matList1)] <- resList1
          matList1[lower.tri(matList1)] <- t(matList1)[lower.tri(matList1)]
          rownames(matList1) <- namsMat
          colnames(matList1) <- namsMat
          resList1 <- matList1
        }
        return(resList1)
      }else {
        return(attr(x[[ionto]][[ilev]][[listNames[1]]][[listNames[2]]]$estimate, "se"))
      }
    })
    names(resLev) <- GOLevel
    return(resLev)
  })
  names(result) <- onto
  return(result)
}

#' Acces to the number of effective bootstrap replicates
#'
#' @export
getNboot <- function(x, ...) {
  UseMethod("getNboot")
}

#' @describeIn getNboot S3 method for class "equivSDhtest"
#' @export
getNboot.equivSDhtest <- function(x) {
  return(attr(x$meth, "nboot"))
}

#' @describeIn getNboot S3 method for class "equivSDhtestList"
#' @export
getNboot.equivSDhtestList <- function(x, simplify = TRUE) {
  result <- lapply(x, function(xi){
    resaux <- lapply(xi, getNboot.equivSDhtest)
    names(resaux) <- names(xi)
    return(resaux)
  })
  namsList1 <- names(result) <- names(x)
  result <- unlist(result)
  if (!simplify) {
    namsMat <- c(names(x[[1]])[1], namsList1)
    nrows <- length(namsList1) + 1
    resMat <- matrix(0.0, nrow = nrows, ncol = nrows)
    resMat[upper.tri(resMat)] <- result
    resMat[lower.tri(resMat)] <- t(resMat)[lower.tri(resMat)]
    rownames(resMat) <- namsMat
    colnames(resMat) <- namsMat
    return(resMat)
  }else{
    return(result)
  }
}
#'
#' @describeIn getNboot S3 method for class "AllEquivSDhtest"
#' @export
getNboot.AllEquivSDhtest <- function(x, onto, GOLevel, listNames,
                                     simplify = TRUE) {
  if (missing(onto)) {
    onto <- names(x)
  }
  if (missing(GOLevel)) {
    GOLevel <- names(x[[1]])
  }
  allLists <- missing(listNames)
  result <- lapply(onto, function(ionto) {
    resLev <- lapply(GOLevel, function(ilev) {
      if (allLists) {
        namsList1 <- names(x[[ionto]][[ilev]])
        namsMat <- c(names(x[[ionto]][[ilev]][[1]])[1], namsList1)
        resList1 <- sapply(namsList1, function(ilist1) {
          namsList2 <- names(x[[ionto]][[ilev]][[ilist1]])
          resList2 <- sapply(namsList2, function(ilist2) {
            return(attr(x[[ionto]][[ilev]][[ilist1]][[ilist2]]$method, "nboot"))
          })
          names(resList2) <- namsList2
          return(resList2)
        })
        names(resList1) <- namsList1
        resList1 <- unlist(resList1, recursive = FALSE)
        if (!simplify) {
          nrows <- length(namsList1) + 1
          matList1 <- matrix(0.0, nrow = nrows, ncol = nrows)
          matList1[upper.tri(matList1)] <- resList1
          matList1[lower.tri(matList1)] <- t(matList1)[lower.tri(matList1)]
          rownames(matList1) <- namsMat
          colnames(matList1) <- namsMat
          resList1 <- matList1
        }
        return(resList1)
      }else {
        return(attr(x[[ionto]][[ilev]][[listNames[1]]][[listNames[2]]]$method, "nboot"))
      }
    })
    names(resLev) <- GOLevel
    # if (simplify) {
    #   resLev <- unlist(resLev, recursive = FALSE)
    # }
    return(resLev)
  })
  names(result) <- onto
  # if (simplify) {
  #   result <- unlist(result, recursive = FALSE)
  # }
  return(result)
}
