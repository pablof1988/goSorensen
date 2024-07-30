#' Studentized Sorensen-Dice dissimilarity statistic
#'
#' Efficient computation of the studentized statistic (^dis - dis) / ^se where 'dis' stands
#' for the "population" value of the Sorensen-Dice dissimilarity, '^dis' for its estimated
#' value and '^se'for the estimate of the standard error of '^dis'. Internally used in
#' bootstrap computations.
#'
#' @param xBoot either an object of class "table", "matrix" or "numeric" representing
#' a 2x2 contingency table of joint enrichment.
#' @param dis the "known" value of the population dissimilarity.
#'
#' @return A numeric value, the result of computing (^dis - dis) / ^se.
#'
#' @details This function is repeatedly evaluated during bootstrap iterations.
#' Given a contingency table 'x' of mutual enrichment (the "true" dataset):
#'
#'\tabular{rr}{
#' \eqn{n_{11}} \tab \eqn{n_{10}} \cr
#' \eqn{n_{01}} \tab \eqn{n_{00}},
#'}
#'
#' summarizing the status of mutual presence of enrichment in two gene lists, where
#' the subindex '11' corresponds to those GO terms enriched in both lists,
#' '01' to terms enriched in the second list but not in the first one,
#' '10' to terms enriched in the first list but not enriched in the second one and
#' '00' to those GO terms non enriched in both gene lists, i.e., to the double negatives.
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
boot.tStat <- function(xBoot, dis) {
  nu <- sum(xBoot[seq_len(3)])
  dBoot <- (xBoot[2] + xBoot[3]) / (nu + xBoot[1])
  p11 <- xBoot[1] / nu
  p11plus <- 1 + p11
  se <- 2 * sqrt(p11 * (1 - p11) / (nu - 1)) / (p11plus * p11plus)
  return((dBoot - dis) / se)
}
