mdes.ird1r1 <- function(power = .80, alpha = .05, two.tailed = TRUE,
                        rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                        r2 = 0, g = 0, p = .50, n) {
  df <- n - g - 2
  d <- .d(p, k1, k2, dists, rhots)
  sse <- sqrt(d * (1 - r2) / (p * (1 - p) * n))
  parms <- as.list(environment())
  .error.handler(parms)
  mdes <- .mdes(power, alpha, sse, df, two.tailed)
  colnames(mdes) <- c("mdes", paste0(100 * (1 - round(alpha, 2)), "% lcl"),
                      paste0(100 * (1 - round(alpha, 2)), "% ucl"))
  .summarize.mdes(power, alpha, sse, df, two.tailed, mdes)
  mdes.out <- list(parms = parms, mdes = mdes)
  class(mdes.out) <- c("mdes", "ird1r1")
  return(invisible(mdes.out))
}

# example
# mdes.ird1r1(n = 400)

power.ird1r1 <- function(es = .25, alpha = .05, two.tailed = TRUE,
                         rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                         r2 = 0, g = 0,
                         p = .50, n) {
  df <- n - g - 2
  d <- .d(p, k1, k2, dists, rhots)
  sse <- sqrt(d * (1 - r2) / (p * (1 - p) * n))
  parms <- as.list(environment())
  .error.handler(parms)
  power <- .power(es, alpha, sse, df, two.tailed)
  .summarize.power(es, alpha, sse, df, two.tailed, power)
  power.out <-  list(parms = parms, power = power)
  class(power.out) <- c("power", "ird1r1")
  return(invisible(power.out))
}
# example
# power.ird1r1(es = 0.466, n = 400)
