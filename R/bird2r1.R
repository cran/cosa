mdes.bird2r1 <- function(power = .80, alpha = .05, two.tailed = TRUE,
                         rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                         rho2, omega2, r21 = 0, r2t2 = 0, g2 = 0,
                         p = .50, n1, n2) {
  df <- n2 - g2 - 1
  d <- .d(p, k1, k2, dists, rhots)
  sse <- sqrt(rho2 * omega2 * (1 - r2t2) / n2 +
                d * (1 - rho2) * (1 - r21) / (p * (1 - p) * n2 * n1))
  parms <- as.list(environment())
  .error.handler(parms)
  mdes <- .mdes(power, alpha, sse, df, two.tailed)
  colnames(mdes) <- c("mdes", paste0(100 * (1 - round(alpha, 2)), "% lcl"),
                      paste0(100 * (1 - round(alpha, 2)), "% ucl"))
  .summarize.mdes(power, alpha, sse, df, two.tailed, mdes)
  mdes.out <- list(parms = parms, mdes = mdes)
  class(mdes.out) <- c("mdes", "bird2r1")
  return(invisible(mdes.out))
}

# example
# mdes.bird2r1(rho2 = .35, omega2 = .10, n1 = 83, n2 = 480)

power.bird2r1 <- function(es = .25, alpha = .05, two.tailed = TRUE,
                          rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                          rho2, omega2, r21 = 0, r2t2 = 0, g2 = 0,
                          p = .50, n1, n2) {
  df <- n2 - g2 - 1
  d <- .d(p, k1, k2, dists, rhots)
  sse <- sqrt(rho2 * omega2 * (1 - r2t2) / n2 +
                d * (1 - rho2) * (1 - r21) / (p * (1 - p) * n2 * n1))
  parms <- as.list(environment())
  .error.handler(parms)
  power <- .power(es, alpha, sse, df, two.tailed)
  .summarize.power(es, alpha, sse, df, two.tailed, power)
  power.out <-  list(parms = parms, power = power)
  class(power.out) <- c("power", "bird2r1")
  return(invisible(power.out))
}
# example
# power.bird2r1(es = 0.0446, rho2 = .35, omega2 = .10, n1 = 83, n2 = 480)

cosa.bird2r1 <- function(cn1 = 0, cn2 = 0, cost = NULL,
                        n1 = NULL, n2 = NULL, p = NULL, n0 = c(10, 100), p0 = .50,
                        constrain = "power", local.solver = c("LBFGS", "SLSQP", "MMA", "COBYLA"),
                        rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                        power = .80, es = .25, alpha = .05, two.tailed = TRUE,
                        rho2, omega2, g2 = 0, r21 = 0, r2t2 = 0) {
  parms <- as.list(environment())
  .error.handler(parms)
  fun <- "cosa.bird2r1"
  lb <- c(1, 1, g2 + 2)

  .df <- quote(n2 - g2 - 1)
  .sse <- quote(sqrt(rho2 * omega2 * (1 - r2t2) / n2 +
                       d * (1 - rho2) * (1 - r21) / (p * (1 - p) * n2 * n1)))
  .cost <- quote(n2 * cn2 +
                   n2 * n1 * (cn1[2] + p * (cn1[1] - cn1[2])))

  cosa.out <- do.call(".cosa", parms)
  class(cosa.out) <- c("cosa", "bird2r1")
  return(invisible(cosa.out))
}
# examples
# unconstrained
# cosa.bird2r1(constrain = "power", rho2 = .20, omega2 = .20, local.solver = "lbfgs")
# constrained
# cosa.bird2r1(cn1 = c(10,5), cn2 = 20, cost = 5000,
#              constrain = "cost", rho2 = .20, omega2 = .20)
