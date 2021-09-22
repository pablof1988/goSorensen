boot.tStat <- function(xBoot, dis) {
  nu <- sum(xBoot[1:3])
  dBoot <- (xBoot[2] + xBoot[3]) / (nu + xBoot[1])
  p11 <- xBoot[1] / nu
  p11plus <- 1 + p11
  se <- 2 * sqrt(p11 * (1 - p11) / (nu - 1)) / (p11plus * p11plus)
  return((dBoot - dis) / se)
}
