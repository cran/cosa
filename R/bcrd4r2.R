mdes.bcrd4r2 <- function(power = .80, alpha = .05, two.tailed = TRUE,
                         rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                         rho2, rho3, rho4, omega3, omega4,
                         r21 = 0, r22 = 0, r2t3 = 0, r2t4 = 0, g4 = 0,
                         p = .50, n1, n2, n3, n4) {

  ss <- c(n1, n2, n3, n4, g4)
  if(any(ss < 0) || !is.numeric(ss) || length(ss) > 5){
    stop("Incorrect value for sample size or number of covariates", call. = FALSE)
  }

  var <- c(r21, r22, r2t3, r2t4, rho2, rho3, rho4, omega3, omega4)
  if(any(var < 0) || any(var > 1) || !is.numeric(var) || length(var) > 9){
    stop("Incorrect value for [0, 1] bounded arguments", call. = FALSE)
  }

  rate <- c(p, power, alpha)
  if(any(rate < .01) || any(rate > .99) || !is.numeric(rate) || length(rate) > 3){
    stop("Incorrect value for [.01, .99] bounded arguments", call. = FALSE)
  }

  if(!is.logical(two.tailed) || length(two.tailed) > 1){
    stop("Non-logical value for 'two.tailed'", call. = FALSE)
  }

  if(any(n4 - g4 < 2)){
    stop("Insufficient sample size, increase 'n4'", call. = FALSE)
  }

  df <- n4 - g4 - 1
  d <- .d(p, k1, k2, dists, rhots)
  sse <- sqrt(rho4 * omega4 * (1 - r2t4) / n4 +
                rho3 * omega3 * (1 - r2t3) / (n4 * n3) +
                d * (rho2 * (1 - r22) / (p * (1 - p) * n4 * n3 * n2) +
                       (1 - rho4 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n4 * n3 * n2 * n1)))

  mdes <- .mdes(power, alpha, sse, df, two.tailed)
  colnames(mdes) <- c("mdes", paste0(100 * (1 - round(alpha, 2)), "%lcl"),
                      paste0(100 * (1 - round(alpha, 2)), "%ucl"))
  mdes.out <- list(parms = list(power = power, alpha = alpha, two.tailed = two.tailed,
                                rhots = rhots, k1 = k1, k2 = k2, dists = dists,
                                rho2 = rho2, rho3 = rho3, rho4 = rho4,
                                omega3 = omega3, omega4 = omega4,
                                r21 = r21, r22 = r22, r2t3 = r2t3, r2t4 = r2t4,
                                g4 = g4, p = p, n1 = n1, n2 = n2, n3 = n3, n4 = n4),
                   df = df,
                   sse = sse,
                   mdes = mdes)
  print(round(mdes, 3))
  class(mdes.out) <- c("mdes", "bcrd4r2")
  return(invisible(mdes.out))
}

# example
# mdes.bcrd4r2(rho4 = .05, rho3 = .15, rho2 = .15, omega4 = .50, omega3 = .50, n1 = 10, n2 = 4, n3 = 4, n4 = 20)

power.bcrd4r2 <- function(es = .25, alpha = .05, two.tailed = TRUE,
                          rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                          rho2, rho3, rho4, omega3, omega4,
                          r21 = 0, r22 = 0, r2t3 = 0, r2t4 = 0, g4 = 0,
                          p = .50, n1, n2, n3, n4) {

  ss <- c(n1, n2, n3, n4, g4)
  if(any(ss < 0) || !is.numeric(ss) || length(ss) > 5){
    stop("Incorrect value for sample size or number of covariates", call. = FALSE)
  }

  var <- c(r21, r22, r2t3, r2t4, rho2, rho3, rho4, omega3, omega4)
  if(any(var < 0) || any(var > 1) || !is.numeric(var) || length(var) > 9){
    stop("Incorrect value for [0, 1] bounded arguments", call. = FALSE)
  }

  rate <- c(p, alpha)
  if(any(rate < .01) || any(rate > .99) || !is.numeric(rate) || length(rate) > 2){
    stop("Incorrect value for [.01, .99] bounded arguments", call. = FALSE)
  }

  if(!is.logical(two.tailed) || length(two.tailed) > 1){
    stop("Non-logical value for 'two.tailed'", call. = FALSE)
  }

  if(any(n4 - g4 < 2)){
    stop("Insufficient sample size, increase 'n4'", call. = FALSE)
  }

  if(any(es <= 0) || !is.numeric(es) || length(es) > 1){
    stop("Incorrect value for 'es'", call. = FALSE)
  }

  df <- n4 - g4 - 1
  d <- .d(p, k1, k2, dists, rhots)
  sse <- sqrt(rho4 * omega4 * (1 - r2t4) / n4 +
                rho3 * omega3 * (1 - r2t3) / (n4 * n3) +
                d * (rho2 * (1 - r22) / (p * (1 - p) * n4 * n3 * n2) +
                       (1 - rho4 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n4 * n3 * n2 * n1)))

  power <- .power(es, alpha, sse, df, two.tailed)
  power.out <-  list(parms = list(es = es, alpha = alpha, two.tailed = two.tailed,
                                  rhots = rhots, k1 = k1, k2 = k2, dists = dists,
                                  rho2 = rho2, rho3 = rho3, rho4 = rho4,
                                  omega3 = omega3, omega4 = omega4,
                                  r21 = r21, r22 = r22, r2t3 = r2t3, r2t4 = r2t4,
                                  g4 = g4, p = p, n1 = n1, n2 = n2, n3 = n3, n4 = n4),
                     df = df,
                     sse = sse,
                     power = power)
  names(power) <- "power"
  print(round(power, 3))
  class(power.out) <- c("power", "bcrd4r2")
  return(invisible(power.out))
}
# example
# power.bcrd4r2(es = .292, rho4 = .05, rho3 = .15, rho2 = .15, omega4 = .50, omega3 = .50, n1 = 10, n2 = 4, n3 = 4, n4 = 20)

cosa.bcrd4r2 <- function(cn1 = 0, cn2 = 0, cn3 = 0, cn4 = 0, cost = NULL,
                         n1 = NULL, n2 = NULL, n3 = NULL, n4 = NULL, p = NULL,
                         n0 = c(10, 3, 100, 5 + g4), p0 = .50,
                         constrain = "power", round = TRUE,
                         local.solver = c("LBFGS", "SLSQP", "MMA", "COBYLA"),
                         rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                         power = .80, es = .25, alpha = .05, two.tailed = TRUE,
                         rho2, rho3, rho4, omega3, omega4,
                         g4 = 0, r21 = 0, r22 = 0, r2t3 = 0, r2t4 = 0) {

  ss <- c(n1, n2, n3, n4, g4)
  if(any(ss < 0) || !is.numeric(ss) || length(ss) > 9){
    stop("Incorrect value for sample size or number of covariates", call. = FALSE)
  }

  var <- c(r21, r22, r2t3, r2t4, rho2, rho3, rho4, omega3, omega4)
  if(any(var < 0) || any(var > 1) || !is.numeric(var) || length(var) > 9){
    stop("Incorrect value for [0, 1] bounded arguments", call. = FALSE)
  }

  rate <- c(alpha, power)
  if(any(rate < .01) || any(rate > .99) || !is.numeric(rate) || length(rate) > 2){
    stop("Incorrect value for [.01, .99] bounded arguments", call. = FALSE)
  }

  if(!is.logical(two.tailed) || length(two.tailed) > 1){
    stop("Non-logical value for 'two.tailed'", call. = FALSE)
  }

  if(any(n4 - g4 < 2)){
    stop("Insufficient sample size, increase 'n4'", call. = FALSE)
  }

  if(any(es <= 0) || !is.numeric(es) || length(es) > 1){
    stop("Incorrect value for 'es'", call. = FALSE)
  }


  fun <- "cosa.bcrd4r2"
  lb <- c(1, 1, 1, g4 + 2)

  .df <- quote(n4 - g4 - 1)
  .sse <- quote(sqrt(rho4 * omega4 * (1 - r2t4) / n4 +
                       rho3 * omega3 * (1 - r2t3) / (n4 * n3) +
                       d * (rho2 * (1 - r22) / (p * (1 - p) * n4 * n3 * n2) +
                              (1 - rho4 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n4 * n3 * n2 * n1))))
  .cost <- quote(n4 * cn4 +
                   n4 * n3 * cn3 +
                   n4 * n3 * n2 * (cn2[2] + p * (cn2[1] - cn2[2])) +
                   n4 * n3 * n2 * n1 * (cn1[2] + p * (cn1[1] - cn1[2])))

  cosa <- .cosa(cn1 = cn1, cn2 = cn2, cn3 = cn3, cn4 = cn4, cost = cost,
                constrain = constrain, round = round, local.solver = local.solver,
                power = power, es = es, alpha = alpha, two.tailed = two.tailed,
                rhots = rhots, k1 = k1, k2 = k2, dists = dists,
                rho2 = rho2, rho3 = rho3, rho4 = rho4,
                omega3 = omega3, omega4 = omega4,
                r21 = r21, r22 = r22, r2t3 = r2t3, r2t4 = r2t4,
                g4 = g4, p0 = p0, p = p, n0 = n0,
                n1 = n1, n2 = n2, n3 = n3, n4 = n4)
  cosa.out <- list(parms = list(cn1 = cn1, cn2 = cn2, cn3 = cn3, cn4 = cn4, cost = cost,
                                constrain = constrain, round = round, local.solver = local.solver,
                                power = power, es = es, alpha = alpha, two.tailed = two.tailed,
                                rhots = rhots, k1 = k1, k2 = k2, dists = dists,
                                rho2 = rho2, rho3 = rho3, rho4 = rho4,
                                omega3 = omega3, omega4 = omega4,
                                r21 = r21, r22 = r22, r2t3 = r2t3, r2t4 = r2t4,
                                g4 = g4, p0 = p0, p = p, n0 = n0,
                                n1 = n1, n2 = n2, n3 = n3, n4 = n4),
                   cosa = cosa)
  print(round(cosa, 3))
  class(cosa.out) <- c("cosa", "bcrd4r2")
  return(invisible(cosa.out))
}
# examples
# unconstrained
# cosa.bcrd4r2(constrain = "power", rho2 = .20, rho3 = .10, rho4 = .05, omega3 = .20, omega4 = .20)
# constrained
# cosa.bcrd4r2(cn1 = c(10,5), cn2 = c(20,10), cn3 = 50, cn4 = 100, cost = 20000,
#              constrain = "cost", rho2 = .20, rho3 = .10, rho4 = .05, omega3 = .20, omega4 = .20)
