mdes.bird4r1 <- function(power = .80, alpha = .05, two.tailed = TRUE,
                         rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                         rho2, rho3, rho4, omega2, omega3, omega4,
                         r21 = 0, r2t2 = 0, r2t3 = 0, r2t4 = 0, g4 = 0,
                         p = .50, n1, n2, n3, n4) {
  df <- n4 - g4 - 1
  d <- .d(p, k1, k2, dists, rhots)
  sse <- sqrt(rho4 * omega4 * (1 - r2t4) / n4 +
                rho3 * omega3 * (1 - r2t3) / (n4 * n3) +
                rho2 * omega2 * (1 - r2t2) / (n4 * n3 * n2) +
                d * (1 - rho4 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n4 * n3 * n2 * n1))
  parms <- as.list(environment())
  .error.handler(parms)
  mdes <- .mdes(power, alpha, sse, df, two.tailed)
  colnames(mdes) <- c("mdes", paste0(100 * (1 - round(alpha, 2)), "% lcl"),
                      paste0(100 * (1 - round(alpha, 2)), "% ucl"))
  .summarize.mdes(power, alpha, sse, df, two.tailed, mdes)
  mdes.out <- list(parms = parms, mdes = mdes)
  class(mdes.out) <- c("mdes", "bird4r1")
  return(invisible(mdes.out))
}
# example
# mdes.bird4r1(power = .80, rho4 = .05, rho3 = .15, rho2 = .15,
#              omega4 = .50, omega3 = .50, omega2 = .50,
#              n1 = 10, n2 = 4, n3 = 27, n4 = 10)

power.bird4r1 <- function(es = .25, alpha = .05, two.tailed = TRUE,
                          rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                          rho2, rho3, rho4, omega2, omega3, omega4,
                          r21 = 0, r2t2 = 0, r2t3 = 0, r2t4 = 0, g4 = 0,
                          p = .50, n1, n2, n3, n4){
  df <- n4 - g4 - 1
  d <- .d(p, k1, k2, dists, rhots)
  sse <- sqrt(rho4 * omega4 * (1 - r2t4) / n4 +
                rho3 * omega3 * (1 - r2t3) / (n4 * n3) +
                rho2 * omega2 * (1 - r2t2) / (n4 * n3 * n2) +
                d * (1 - rho4 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n4 * n3 * n2 * n1))
  parms <- as.list(environment())
  .error.handler(parms)
  power <- .power(es, alpha, sse, df, two.tailed)
  .summarize.power(es, alpha, sse, df, two.tailed, power)
  power.out <-  list(parms = parms, power = power)
  class(power.out) <- c("power", "bird4r1")
  return(invisible(power.out))
}
# example
# power.bird4r1(es = 0.1863523, rho4 = .05, rho3 = .15, rho2 = .15,
#              omega4 = .50, omega3 = .50, omega2 = .50,
#              n1 = 10, n2 = 4, n3 = 27, n4 = 10)

cosa.bird4r1 <- function(cn1 = 0, cn2 = 0, cn3 = 0, cn4 = 0, cost = NULL,
                         n1 = NULL, n2 = NULL, n3 = NULL, n4 = NULL, p = NULL,
                         n0 = c(10, 3, 100, 10), p0 = .50,
                         constrain = "power", local.solver = c("LBFGS", "SLSQP", "MMA", "COBYLA"),
                         rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                         power = .80, es = .25, alpha = .05, two.tailed = TRUE,
                         rho2, rho3, rho4, omega2, omega3, omega4,
                         g4 = 0, r21 = 0, r2t2 = 0, r2t3 = 0, r2t4 = 0) {
  parms <- as.list(environment())
  .error.handler(parms)
  fun <- "cosa.bird4r1"
  lb <- c(1, 1, 1, g4 + 2)

  .df <- quote(n4 - g4 - 1)
  .sse <- quote(sqrt(rho4 * omega4 * (1 - r2t4) / n4 +
                       rho3 * omega3 * (1 - r2t3) / (n4 * n3) +
                       rho2 * omega2 * (1 - r2t2) / (n4 * n3 * n2) +
                       d * (1 - rho4 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n4 * n3 * n2 * n1)))
  .cost <- quote(n4 * cn4 +
                   n4 * n3 * cn3 +
                   n4 * n3 * n2 * cn2 +
                   n4 * n3 * n2 * n1 * (cn1[2] + p * (cn1[1] - cn1[2])))

  cosa.out <- do.call(".cosa", parms)
  class(cosa.out) <- c("cosa", "bird4r1")
  return(invisible(cosa.out))
}
# examples
# unconstrained
# cosa.bird4r1(constrain = "power", rho2 = .20, rho3 = .10, rho4 = .05, omega2 = .20, omega3 = .20, omega4 = .20)
# constrained
# cosa.bird4r1(cn1 = c(10,5), cn2 = 20, cn3 = 50, cn4 = 100, cost = 20000,
#              constrain = "cost", rho2 = .20, rho3 = .10, rho4 = .05, omega2 = .20, omega3 = .20, omega4 = .20)
