#' Studentized Sorensen-Dice dissimilarity statistic
#' 
#' Efficient computation of the studentized statistic (^dis - dis) / ^se where 'dis' stands
#' for the "population" value of the Sorensen-Dice dissimilarity, '^dis' for its estimated
#' value and '^se'for the estimate of the standard error of '^dis'. Internally used in
#' bootstrap computations.
#'
#' @param xBoot either an object of class "table", "matrix" or "numeric" representing
#' a 2x2 contingency table of joint enrichment.
#' @param dis the "knwon" value of the population dissimilarity.
#'
#' @return A numeric value, the result of computing (^dis - dis) / ^se.
#'
#' @details This function is repeatedly evaluated during bootstrap iterations.
#' Given a contingency table 'x' of mutual enrichment (the "true" dataset):
#'
#'\tabular{rr}{
#' n_11 \tab n_01 \cr
#' n_10 \tab n_00,
#'}
#'
#' summarizing the status of mutual presence of enrichment in two gene lists, where
#' the subindex '11' corresponds to those GO items enriched in both lists,
#' '01' to items enriched in the second list but not in the first one,
#' '10' to items enriched in the first list but not enriched in the second one and
#' '00' to those GO items non enriched in both gene lists, i.e., to the double negatives.
#'
#' A typical bootstrap iteration consists in repeatedly generating four frequencies
#' from a multinomial of parameters size = sum(n_ij), i,j = 1, 0 and probabilities
#' (n_11/size, n_10/size, n_10/size, n_00/size).
#' The argument 'xBoot' corresponds to each one of these bootstrap resamples (indiferenly
#' represented in form of a 2x2 "table" or "matrix" or as a numeric vector)
#' In each bootstrap iteration, the value of the "true" known 'dis' is the dissimilarity
#' which was computed from 'x' (a constant, known value in the full iteration) and the
#' values of '^dis' and '^se' are internally computed from the bootstrap data 'xBoot'.
#'
#' @examples
#' # Gene lists 'atlas' and 'sanger' in 'Cangenes' dataset. Table of joint enrichment
#' # of GO items in ontology BP at level 3.
#' tab_atlas.sanger_BP3
#' # Generate B = 100 bootstrap values of the studentized statistic:
#' # 1) Sorensen-Dice dissimilarity for the "real" data 'tab_atlas.sanger_BP3'.
#' d <- dSorensen(tab_atlas.sanger_BP3)
#' # 2) 100 multinomial data resamples:
#' size <- sum(tab_atlas.sanger_BP3)
#' set.seed(12345)
#' xBoots <- rmultinom(100, size, tab_atlas.sanger_BP3 / size)
#' # Each colum of 'xBoots' corresponds to a frequencies boostrap resample:
#' # 3) 100 bootstrap studentized statistics:
#' bootVals <- apply(xBoots, 2, boot.tStat, dis = d)
#' bootVals

boot.tStat <- function(xBoot, dis) {
  nu <- sum(xBoot[1:3])
  dBoot <- (xBoot[2] + xBoot[3]) / (nu + xBoot[1])
  p11 <- xBoot[1] / nu
  p11plus <- 1 + p11
  se <- 2 * sqrt(p11 * (1 - p11) / (nu - 1)) / (p11plus * p11plus)
  return((dBoot - dis) / se)
}
