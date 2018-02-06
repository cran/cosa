mdes.crd3r3 <- function(power = .80, alpha = .05, two.tailed = TRUE,
                        rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                        rho2, rho3, r21 = 0, r22 = 0, r23 = 0, g3 = 0,
                        p = .50, n1, n2, n3) {
  df <- n3 - g3 - 2
  d <- .d(p, k1, k2, dists, rhots)
  sse<- sqrt(d * (rho3 * (1 - r23) / (p * (1 - p) * n3) +
                rho2 * (1 - r22) / (p * (1 - p) * n3 * n2) +
               (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n3 * n2 * n1)))
  parms <- as.list(environment())
  .error.handler(parms)
  mdes <- .mdes(power, alpha, sse, df, two.tailed)
  colnames(mdes) <- c("mdes", paste0(100 * (1 - round(alpha, 2)), "% lcl"),
                      paste0(100 * (1 - round(alpha, 2)), "% ucl"))
  .summarize.mdes(power, alpha, sse, df, two.tailed, mdes)
  mdes.out <- list(parms = parms, mdes = mdes)
  class(mdes.out) <- c("mdes", "crd3r3")
  return(invisible(mdes.out))
}

# example
# mdes.crd3r3(rho3 = .06, rho2 = .17, n1 = 15, n2 = 3, n3 = 60)


power.crd3r3 <- function(es = .25, alpha = .05, two.tailed = TRUE,
                         rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                         rho2, rho3, r21 = 0, r22 = 0, r23 = 0, g3 = 0,
                         p = .50, n1, n2, n3) {
  df <- n3 - g3 - 2
  d <- .d(p, k1, k2, dists, rhots)
  sse <- sqrt(d * (rho3 * (1 - r23) / (p * (1 - p) * n3) +
                    rho2 * (1 - r22) / (p * (1 - p) * n3 * n2) +
                    (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n3 * n2 * n1)))
  parms <- as.list(environment())
  .error.handler(parms)
  power <- .power(es, alpha, sse, df, two.tailed)
  .summarize.power(es, alpha, sse, df, two.tailed, power)
  power.out <-  list(parms = parms, power = power)
  class(power.out) <- c("power", "crd3r3")
  return(invisible(power.out))
}
# example
# power.crd3r3(es = .446, rho3 = .06, rho2 = .17, n1 = 15, n2 = 3, n3 = 60)

cosa.crd3r3 <- function(cn1 = 0, cn2 = 0, cn3 = 0, cost = NULL,
                        n1 = NULL, n2 = NULL, n3 = NULL, p = NULL, n0 = c(10, 3, 100), p0 = .50,
                        constrain = "power", local.solver = c("LBFGS", "SLSQP", "MMA", "COBYLA"),
                        rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                        power = .80, es = .25, alpha = .05, two.tailed = TRUE,
                        rho2, rho3, g3 = 0, r21 = 0, r22 = 0, r23 = 0) {
  parms <- as.list(environment())
  .error.handler(parms)
  fun <- "cosa.crd3r3"
  lb <- c(1, 1, g3 + 3)

  .df <- quote(n3 - g3 - 2)
  .sse <- quote(sqrt(d * (rho3 * (1 - r23) / (p * (1 - p) * n3) +
                            rho2 * (1 - r22) / (p * (1 - p) * n3 * n2) +
                            (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n3 * n2 * n1))))
  .cost <- quote(n3 * (cn3[2] + p * (cn3[1] - cn3[2])) +
                   n3 * n2 * (cn2[2] + p * (cn2[1] - cn2[2])) +
                   n3 * n2 * n1 * (cn1[2] + p * (cn1[1] - cn1[2])))

  cosa.out <- do.call(".cosa", parms)
  class(cosa.out) <- c("cosa", "crd3r3")
  return(invisible(cosa.out))
}
# examples
# unconstrained
# cosa.crd3r3(constrain = "power", rho2 = .20, rho3 = .10)
# constrained
# cosa.crd3r3(cn1 = c(10,5), cn2 = c(20,10), cn3 = c(50,25), cost = 10000,
#             constrain = "cost", rho2 = .20, rho3 = .10)
