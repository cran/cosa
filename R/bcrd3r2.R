mdes.bcrd3r2 <- function(power = .80, alpha = .05, two.tailed = TRUE,
                         rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                         rho2, rho3, omega3, r21 = 0, r22 = 0, r2t3 = 0, g3 = 0,
                         p = .50, n1, n2, n3) {
  df <- n3 - g3 - 1
  d <- .d(p, k1, k2, dists, rhots)
  sse <- sqrt(rho3 * omega3 * (1 - r2t3) / n3 +
                d * (rho2 * (1 - r22) / (p * (1 - p) * n2 * n3) +
                     (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n3 * n2 * n1)))
  parms <- as.list(environment())
  .error.handler(parms)
  mdes <- .mdes(power, alpha, sse, df, two.tailed)
  colnames(mdes) <- c("mdes", paste0(100 * (1 - round(alpha, 2)), "% lcl"),
                      paste0(100 * (1 - round(alpha, 2)), "% ucl"))
  .summarize.mdes(power, alpha, sse, df, two.tailed, mdes)
  mdes.out <- list(parms = parms, mdes = mdes)
  class(mdes.out) <- c("mdes", "bcrd3r2")
  return(invisible(mdes.out))
}

# example
# mdes.bcrd3r2(rho2 = .10, rho3 = .20, omega3 = .30, n1 = 20, n2 = 44, n3 = 10)


power.bcrd3r2 <- function(es = .25, alpha = .05, two.tailed = TRUE,
                          rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                          rho2, rho3, omega3, r21 = 0, r22 = 0, r2t3 = 0, g3 = 0,
                          p = .50, n1, n2, n3) {
  df <- n3 - g3 - 1
  d <- .d(p, k1, k2, dists, rhots)
  sse <- sqrt(rho3 * omega3 * (1 - r2t3) / n3 +
                d * (rho2 * (1 - r22) / (p * (1 - p) * n2 * n3) +
                       (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n3 * n2 * n1)))
  parms <- as.list(environment())
  .error.handler(parms)
  power <- .power(es, alpha, sse, df, two.tailed)
  .summarize.power(es, alpha, sse, df, two.tailed, power)
  power.out <-  list(parms = parms, power = power)
  class(power.out) <- c("power", "bcrd3r2")
  return(invisible(power.out))
}
# example
# power.bcrd3r2(es = .305, rho2 = .10, rho3 = .20, omega3 = .30, n1 = 20, n2 = 44, n3 = 10)

cosa.bcrd3r2 <- function(cn1 = 0, cn2 = 0, cn3 = 0, cost = NULL,
                        n1 = NULL, n2 = NULL, n3 = NULL, p = NULL, n0 = c(10, 3, 100), p0 = .50,
                        constrain = "power", local.solver = c("LBFGS", "SLSQP", "MMA", "COBYLA"),
                        rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                        power = .80, es = .25, alpha = .05, two.tailed = TRUE,
                        rho2, rho3, omega3, g3 = 0, r21 = 0, r22 = 0, r2t3 = 0) {
  parms <- as.list(environment())
  .error.handler(parms)
  fun <- "cosa.bcrd3r2"
  lb <- c(1, 1, g3 + 2)

  .df <- quote(n3 - g3 - 1)
  .sse <- quote(sqrt(rho3 * omega3 * (1 - r2t3) / n3 +
                       d * (rho2 * (1 - r22) / (p * (1 - p) * n2 * n3) +
                              (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n3 * n2 * n1))))
  .cost <- quote(n3 * cn3 +
                   n3 * n2 * (cn2[2] + p * (cn2[1] - cn2[2])) +
                   n3 * n2 * n1 * (cn1[2] + p * (cn1[1] - cn1[2])))

  cosa.out <- do.call(".cosa", parms)
  class(cosa.out) <- c("cosa", "bcrd3r2")
  return(invisible(cosa.out))
}
# examples
# unconstrained
# cosa.bcrd3r2(constrain = "power", rho2 = .20, rho3 = .10,  omega3 = .20)
# constrained
# cosa.bcrd3r2(cn1 = c(10,5), cn2 = c(20,10), cn3 = 50, cost = 10000,
#              constrain = "cost", rho2 = .20, rho3 = .10, omega3 = .20)
