mdes.crd4r4 <- function(power = .80, alpha = .05, two.tailed = TRUE,
                        rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                        rho2, rho3, rho4, r21 = 0, r22 = 0, r23 = 0, r24 = 0, g4 = 0,
                        p = .50, n1, n2, n3, n4) {
  df <- n4 - g4 - 2
  d <- .d(p, k1, k2, dists, rhots)
  sse <- sqrt(d * (rho4 * (1 - r24) / (p * (1 - p) * n4) +
                     rho3 * (1 - r23) / (p * (1 - p) * n4 * n3) +
                     rho2 * (1 - r22) / (p * (1 - p) * n4 * n3 * n2) +
                     (1 - rho4 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n4 * n3 * n2 * n1)))
  parms <- as.list(environment())
  .error.handler(parms)
  mdes <- .mdes(power, alpha, sse, df, two.tailed)
  colnames(mdes) <- c("mdes", paste0(100 * (1 - round(alpha, 2)), "% lcl"),
                      paste0(100 * (1 - round(alpha, 2)), "% ucl"))
  .summarize.mdes(power, alpha, sse, df, two.tailed, mdes)
  mdes.out <- list(parms = parms, mdes = mdes)
  class(mdes.out) <- c("mdes", "crd4r4")
  return(invisible(mdes.out))
}

# example
# mdes.crd4r4(rho4 = .05, rho3 = .05, rho2 = .10, n1 = 10, n2 = 2, n3 = 3, n4 = 20)

power.crd4r4 <- function(es = .25, alpha = .05, two.tailed = TRUE,
                         rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                         rho2, rho3, rho4, r21 = 0, r22 = 0, r23 = 0, r24 = 0, g4 = 0,
                         p = .50, n1, n2, n3, n4) {
  df <- n4 - g4 - 2
  d <- .d(p, k1, k2, dists, rhots)
  sse <- sqrt(d * (rho4 * (1 - r24) / (p * (1 - p) * n4) +
                     rho3 * (1 - r23) / (p * (1 - p) * n4 * n3) +
                     rho2 * (1 - r22) / (p * (1 - p) * n4 * n3 * n2) +
                     (1 - rho4 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n4 * n3 * n2 * n1)))
  parms <- as.list(environment())
  .error.handler(parms)
  power <- .power(es, alpha, sse, df, two.tailed)
  .summarize.power(es, alpha, sse, df, two.tailed, power)
  power.out <-  list(parms = parms, power = power)
  class(power.out) <- c("power", "crd4r4")
  return(invisible(power.out))
}
# example
# power.crd4r4(es = .683, rho4 = .05, rho3 = .05, rho2 = .10, n1 = 10, n2 = 2, n3 = 3, n4 = 20)

cosa.crd4r4 <- function(cn1 = 0, cn2 = 0, cn3 = 0, cn4 = 0, cost = NULL,
                        n1 = NULL, n2 = NULL, n3 = NULL, n4 = NULL, p = NULL, n0 = c(10, 3, 100, 10), p0 = .50,
                        constrain = "power", local.solver = c("LBFGS", "SLSQP", "MMA", "COBYLA"),
                        rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                        power = .80, es = .25, alpha = .05, two.tailed = TRUE,
                        rho2, rho3, rho4, g4 = 0, r21 = 0, r22 = 0, r23 = 0, r24 = 0) {
  parms <- as.list(environment())
  .error.handler(parms)
  fun <- "cosa.crd4r4"
  lb <- c(1, 1, 1, g4 + 3)

  .df <- quote(n4 - g4 - 2)
  .sse <- quote(sqrt(d * (rho4 * (1 - r24) / (p * (1 - p) * n4) +
                            rho3 * (1 - r23) / (p * (1 - p) * n4 * n3) +
                            rho2 * (1 - r22) / (p * (1 - p) * n4 * n3 * n2) +
                            (1 - rho4 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n4 * n3 * n2 * n1))))
  .cost <- quote(n4 * (cn4[2] + p * (cn4[1] - cn4[2])) +
                   n4 * n3 * (cn3[2] + p * (cn3[1] - cn3[2])) +
                   n4 * n3 * n2 * (cn2[2] + p * (cn2[1] - cn2[2])) +
                   n4 * n3 * n2 * n1 * (cn1[2] + p * (cn1[1] - cn1[2])))

  cosa.out <- do.call(".cosa", parms)
  class(cosa.out) <- c("cosa", "crd4r4")
  return(invisible(cosa.out))
}
# examples
# unconstrained
# cosa.crd4r4(constrain = "power", rho2 = .20, rho3 = .10, rho4 = .05)
# constrained
# cosa.crd4r4(cn1 = c(10,5), cn2 = c(20,10), cn3 = c(50,25), cn4 = c(100,50), cost = 20000,
#             constrain = "cost", rho2 = .20, rho3 = .10, rho4 = .05)
