mdes.bird3r1 <- function(power = .80, alpha = .05, two.tailed = TRUE,
                         rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                         rho2, rho3, omega2, omega3,
                         r21 = 0, r2t2 = 0, r2t3 = 0, g3 = 0,
                         p = .50, n1, n2, n3) {
  df <- n3 - g3 - 1
  d <- .d(p, k1, k2, dists, rhots)
  sse <- sqrt(rho3 * omega3 * (1 - r2t3) / n3 +
                rho2 * omega2 * (1 - r2t2) / (n3 * n2) +
                d * (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n3 * n2 * n1))
  parms <- as.list(environment())
  .error.handler(parms)
  mdes <- .mdes(power, alpha, sse, df, two.tailed)
  colnames(mdes) <- c("mdes", paste0(100 * (1 - round(alpha, 2)), "% lcl"),
                      paste0(100 * (1 - round(alpha, 2)), "% ucl"))
  .summarize.mdes(power, alpha, sse, df, two.tailed, mdes)
  mdes.out <- list(parms = parms, mdes = mdes)
  class(mdes.out) <- c("mdes", "bird3r1")
  return(invisible(mdes.out))
}

# example
# mdes.bird3r1(rho3 = .20, rho2 = .15, omega3 = .10, omega2 = .10, n1 = 69, n2 = 10, n3 = 100)


power.bird3r1 <- function(es = .25, alpha = .05, two.tailed = TRUE,
                          rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                          rho2, rho3, omega2, omega3,
                          r21 = 0, r2t2 = 0, r2t3 = 0, g3 = 0,
                          p = .50, n1, n2, n3) {
  df <- n3 - g3 - 1
  d <- .d(p, k1, k2, dists, rhots)
  sse <- sqrt(rho3 * omega3 * (1 - r2t3) / n3 +
                rho2 * omega2 * (1 - r2t2) / (n3 * n2) +
                d * (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n3 * n2 * n1))
  parms <- as.list(environment())
  .error.handler(parms)
  power <- .power(es, alpha, sse, df, two.tailed)
  .summarize.power(es, alpha, sse, df, two.tailed, power)
  power.out <-  list(parms = parms, power = power)
  class(power.out) <- c("power", "bird3r1")
  return(invisible(power.out))
}
# example
# power.bird3r1(es = 0.05051213, rho3 = .20, rho2 = .15, omega3 = .10, omega2 = .10, n1 = 69, n2 = 10, n3 = 100)

cosa.bird3r1 <- function(cn1 = 0, cn2 = 0, cn3 = 0, cost = NULL,
                         n1 = NULL, n2 = NULL, n3 = NULL, p = NULL,
                         n0 = c(10, 3, 100), p0 = .50,
                         constrain = "power", local.solver = c("LBFGS", "SLSQP", "MMA", "COBYLA"),
                         rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                         power = .80, es = .25, alpha = .05, two.tailed = TRUE,
                         rho2, rho3, omega2, omega3,
                         g3 = 0, r21 = 0, r2t2 = 0, r2t3 = 0) {
  parms <- as.list(environment())
  .error.handler(parms)
  fun <- "cosa.bird3r1"
  lb <- c(1, 1, g3 + 2)

  .df <- quote(n3 - g3 - 1)
  .sse <- quote(sqrt(rho3 * omega3 * (1 - r2t3) / n3 +
                       rho2 * omega2 * (1 - r2t2) / (n3 * n2) +
                       d * (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n3 * n2 * n1)))
  .cost <- quote(n3 * cn3 +
                   n3 * n2 * cn2 +
                   n3 * n2 * n1 * (cn1[2] + p * (cn1[1] - cn1[2])))

  cosa.out <- do.call(".cosa", parms)
  class(cosa.out) <- c("cosa", "bird3r1")
  return(invisible(cosa.out))
}
# examples
# unconstrained
# cosa.bird3r1(constrain = "power", rho2 = .20, rho3 = .10, omega2 = .20, omega3 = .20)
# constrained
# cosa.bird3r1(cn1 = c(10,5), cn2 = 20, cn3 = 50, cost = 10000,
#              constrain = "cost", rho2 = .20, rho3 = .10, omega2 = .20, omega3 = .20)
