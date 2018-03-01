mdes.bird3r1 <- function(power = .80, alpha = .05, two.tailed = TRUE,
                         rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                         rho2, rho3, omega2, omega3,
                         r21 = 0, r2t2 = 0, r2t3 = 0, g3 = 0,
                         p = .50, n1, n2, n3) {

  ss <- c(n1, n2, n3, g3)
  if(any(ss < 0) || !is.numeric(ss) || length(ss) > 4){
    stop("Incorrect value for sample size or number of covariates", call. = FALSE)
  }

  var <- c(r21, r2t2, r2t3, rho2, rho3, omega2, omega3)
  if(any(var < 0) || any(var > 1) || !is.numeric(var) || length(var) > 7){
    stop("Incorrect value for [0, 1] bounded arguments", call. = FALSE)
  }

  rate <- c(p, power, alpha)
  if(any(rate < .01) || any(rate > .99) || !is.numeric(rate) || length(rate) > 3){
    stop("Incorrect value for [.01, .99] bounded arguments", call. = FALSE)
  }

  if(!is.logical(two.tailed) || length(two.tailed) > 1){
    stop("Non-logical value for 'two.tailed'", call. = FALSE)
  }

  if(any(n3 - g3 < 2)){
    stop("Insufficient sample size, increase 'n3'", call. = FALSE)
  }


  df <- n3 - g3 - 1
  d <- .d(p, k1, k2, dists, rhots)
  sse <- sqrt(rho3 * omega3 * (1 - r2t3) / n3 +
                rho2 * omega2 * (1 - r2t2) / (n3 * n2) +
                d * (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n3 * n2 * n1))

  mdes <- .mdes(power, alpha, sse, df, two.tailed)
  colnames(mdes) <- c("mdes", paste0(100 * (1 - round(alpha, 2)), "%lcl"),
                      paste0(100 * (1 - round(alpha, 2)), "%ucl"))
  mdes.out <- list(parms = list(power = power, alpha = alpha, two.tailed = two.tailed,
                                rhots = rhots, k1 = k1, k2 = k2, dists = dists,
                                rho2 = rho2, rho3 = rho3, omega2 = omega2, omega3 = omega3,
                                r21 = r21, r2t2 = r2t2, r2t3 = r2t3,
                                g3 = g3, p = p, n1 = n1, n2 = n2, n3 = n3),
                   df = df,
                   sse = sse,
                   mdes = mdes)
  print(round(mdes, 3))
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

  ss <- c(n1, n2, n3, g3)
  if(any(ss < 0) || !is.numeric(ss) || length(ss) > 4){
    stop("Incorrect value for sample size or number of covariates", call. = FALSE)
  }

  var <- c(r21, r2t2, r2t3, rho2, rho3, omega2, omega3)
  if(any(var < 0) || any(var > 1) || !is.numeric(var) || length(var) > 7){
    stop("Incorrect value for [0, 1] bounded arguments", call. = FALSE)
  }

  rate <- c(p, alpha)
  if(any(rate < .01) || any(rate > .99) || !is.numeric(rate) || length(rate) > 2){
    stop("Incorrect value for [.01, .99] bounded arguments", call. = FALSE)
  }

  if(!is.logical(two.tailed) || length(two.tailed) > 1){
    stop("Non-logical value for 'two.tailed'", call. = FALSE)
  }

  if(any(n3 - g3 < 2)){
    stop("Insufficient sample size, increase 'n3'", call. = FALSE)
  }

  if(any(es <= 0) || !is.numeric(es) || length(es) > 1){
    stop("Incorrect value for 'es'", call. = FALSE)
  }

  df <- n3 - g3 - 1
  d <- .d(p, k1, k2, dists, rhots)
  sse <- sqrt(rho3 * omega3 * (1 - r2t3) / n3 +
                rho2 * omega2 * (1 - r2t2) / (n3 * n2) +
                d * (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n3 * n2 * n1))

  power <- .power(es, alpha, sse, df, two.tailed)
  power.out <-  list(parms = list(es = es, alpha = alpha, two.tailed = two.tailed,
                                  rhots = rhots, k1 = k1, k2 = k2, dists = dists,
                                  rho2 = rho2, rho3 = rho3, omega2 = omega2, omega3 = omega3,
                                  r21 = r21, r2t2 = r2t2, r2t3 = r2t3,
                                  g3 = g3, p = p, n1 = n1, n2 = n2, n3 = n3),
                     df = df,
                     sse = sse,
                     power = power)
  names(power) <- "power"
  print(round(power, 3))
  class(power.out) <- c("power", "bird3r1")
  return(invisible(power.out))
}
# example
# power.bird3r1(es = 0.05051213, rho3 = .20, rho2 = .15, omega3 = .10, omega2 = .10, n1 = 69, n2 = 10, n3 = 100)

cosa.bird3r1 <- function(cn1 = 0, cn2 = 0, cn3 = 0, cost = NULL,
                         n1 = NULL, n2 = NULL, n3 = NULL, p = NULL,
                         n0 = c(10, 3, 100 + g3), p0 = .50,
                         constrain = "power", round = TRUE,
                         local.solver = c("LBFGS", "SLSQP", "MMA", "COBYLA"),
                         rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                         power = .80, es = .25, alpha = .05, two.tailed = TRUE,
                         rho2, rho3, omega2, omega3,
                         g3 = 0, r21 = 0, r2t2 = 0, r2t3 = 0) {

  ss <- c(n1, n2, n3, g3)
  if(any(ss < 0) || !is.numeric(ss) || length(ss) > 7){
    stop("Incorrect value for sample size or number of covariates", call. = FALSE)
  }

  var <- c(r21, r2t2, r2t3, rho2, rho3, omega2, omega3)
  if(any(var < 0) || any(var > 1) || !is.numeric(var) || length(var) > 7){
    stop("Incorrect value for [0, 1] bounded arguments", call. = FALSE)
  }

  rate <- c(alpha, power)
  if(any(rate < .01) || any(rate > .99) || !is.numeric(rate) || length(rate) > 2){
    stop("Incorrect value for [.01, .99] bounded arguments", call. = FALSE)
  }

  if(!is.logical(two.tailed) || length(two.tailed) > 1){
    stop("Non-logical value for 'two.tailed'", call. = FALSE)
  }

  if(any(n3 - g3 < 2)){
    stop("Insufficient sample size, increase 'n3'", call. = FALSE)
  }

  if(any(es <= 0) || !is.numeric(es) || length(es) > 1){
    stop("Incorrect value for 'es'", call. = FALSE)
  }

  fun <- "cosa.bird3r1"
  lb <- c(1, 1, g3 + 2)

  .df <- quote(n3 - g3 - 1)
  .sse <- quote(sqrt(rho3 * omega3 * (1 - r2t3) / n3 +
                       rho2 * omega2 * (1 - r2t2) / (n3 * n2) +
                       d * (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n3 * n2 * n1)))
  .cost <- quote(n3 * cn3 +
                   n3 * n2 * cn2 +
                   n3 * n2 * n1 * (cn1[2] + p * (cn1[1] - cn1[2])))

  cosa <- .cosa(cn1 = cn1, cn2 = cn2, cn3 = cn3, cost = cost,
                constrain = constrain, round = round, local.solver = local.solver,
                power = power, es = es, alpha = alpha, two.tailed = two.tailed,
                rhots = rhots, k1 = k1, k2 = k2, dists = dists,
                rho2 = rho2, rho3 = rho3,
                omega2 = omega2, omega3 = omega3,
                r21 = r21, r2t2 = r2t2, r2t3 = r2t3,
                g3 = g3, p0 = p0, p = p, n0 = n0,
                n1 = n1, n2 = n2, n3 = n3)
  cosa.out <- list(parms = list(cn1 = cn1, cn2 = cn2, cn3 = cn3, cost = cost,
                                constrain = constrain, round = round, local.solver = local.solver,
                                power = power, es = es, alpha = alpha, two.tailed = two.tailed,
                                rhots = rhots, k1 = k1, k2 = k2, dists = dists,
                                rho2 = rho2, rho3 = rho3,
                                omega2 = omega2, omega3 = omega3,
                                r21 = r21, r2t2 = r2t2, r2t3 = r2t3,
                                g3 = g3, p0 = p0, p = p, n0 = n0,
                                n1 = n1, n2 = n2, n3 = n3),
                   cosa = cosa)
  print(round(cosa, 3))
  class(cosa.out) <- c("cosa", "bird3r1")
  return(invisible(cosa.out))
}
# examples
# unconstrained
# cosa.bird3r1(constrain = "power", rho2 = .20, rho3 = .10, omega2 = .20, omega3 = .20)
# constrained
# cosa.bird3r1(cn1 = c(10,5), cn2 = 20, cn3 = 50, cost = 10000,
#              constrain = "cost", rho2 = .20, rho3 = .10, omega2 = .20, omega3 = .20)
