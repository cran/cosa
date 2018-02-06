mdes.crd2r2 <- function(power = .80, alpha = .05, two.tailed = TRUE,
                        rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                        rho2, r21 = 0, r22 = 0, g2 = 0,
                        p = .50, n1, n2) {
  df <- n2 - g2 - 2
  d <- .d(p, k1, k2, dists, rhots)
  sse <- sqrt(d * (rho2 * (1 - r22) / (p * (1 - p) * n2) +
               (1 - rho2) * (1 - r21) / (p * (1 - p) * n2 * n1)))
  parms <- as.list(environment())
  .error.handler(parms)
  mdes <- .mdes(power, alpha, sse, df, two.tailed)
  colnames(mdes) <- c("mdes", paste0(100 * (1 - round(alpha, 2)), "% lcl"),
                      paste0(100 * (1 - round(alpha, 2)), "% ucl"))
  .summarize.mdes(power, alpha, sse, df, two.tailed, mdes)
  mdes.out <- list(parms = parms, mdes = mdes)
  class(mdes.out) <- c("mdes", "crd2r2")
  return(invisible(mdes.out))
}

# example
# mdes.crd2r2(rho2=.20, n1 = 4, n2 = 20)


power.crd2r2 <- function(es = .25, alpha = .05, two.tailed = TRUE,
                         rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                         rho2, r21 = 0, r22 = 0, g2 = 0,
                         p = .50, n1, n2) {
  df <- n2 - g2 - 2
  d <- .d(p, k1, k2, dists, rhots)
  sse <- sqrt(d * (rho2 * (1 - r22) / (p * (1 - p) * n2) +
                     (1 - rho2) * (1 - r21) / (p * (1 - p) * n2 * n1)))
  parms <- as.list(environment())
  .error.handler(parms)
  power <- .power(es, alpha, sse, df, two.tailed)
  .summarize.power(es, alpha, sse, df, two.tailed, power)
  power.out <-  list(parms = parms, power = power)
  class(power.out) <- c("power", "crd2r2")
  return(invisible(power.out))
}
# example
# power.crd2r2(es = 1.391, rho2 = .20, n1 = 4, n2 = 20)

cosa.crd2r2 <- function(cn1 = 0, cn2 = 0, cost = NULL,
                        n1 = NULL, n2 = NULL, p = NULL, n0 = c(10, 100), p0 = .50,
                        constrain = "power", local.solver = c("LBFGS", "SLSQP", "MMA", "COBYLA"),
                        rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                        power = .80, es = .25, alpha = .05, two.tailed = TRUE,
                        rho2, g2 = 0, r21 = 0, r22 = 0) {
  parms <- as.list(environment())
  .error.handler(parms)
  fun <- "cosa.crd2r2"
  lb <- c(1, g2 + 3)

  .df <- quote(n2 - g2 - 2)
  .sse <- quote(sqrt(d * (rho2 * (1 - r22) / (p * (1 - p) * n2) +
                            (1 - rho2) * (1 - r21) / (p * (1 - p) * n2 * n1))))
  .cost <- quote(n2 * (cn2[2] + p * (cn2[1] - cn2[2])) +
                   n2 * n1 * (cn1[2] + p * (cn1[1] - cn1[2])))

  cosa.out <- do.call(".cosa", parms)
  class(cosa.out) <- c("cosa", "crd2r2")
  return(invisible(cosa.out))
}
# examples
# unconstrained
# cosa.crd2r2(constrain = "power", rho2 = .20, local.solver = "lbfgs")
# constrained
# cosa.crd2r2(cn1 = c(10,5), cn2 = c(20,10), cost = 5000,
#             constrain = "cost", rho2 = .20)
