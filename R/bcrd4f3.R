mdes.bcrd4f3 <- function(power = .80, alpha = .05, two.tailed = TRUE,
                        rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                        rho2, rho3, r21 = 0, r22 = 0, r23 = 0, g3 = 0,
                        p = .50, n1, n2, n3, n4) {

  ss <- c(n1, n2, n3, n4, g3)
  if(any(ss < 0) || !is.numeric(ss) || length(ss) > 5){
    stop("Incorrect value for sample size or number of covariates", call. = FALSE)
  }

  var <- c(r21, r22, r23, rho2, rho3)
  if(any(var < 0) || any(var > 1) || !is.numeric(var) || length(var) > 5){
    stop("Incorrect value for [0, 1] bounded arguments", call. = FALSE)
  }

  rate <- c(p, power, alpha)
  if(any(rate < .01) || any(rate > .99) || !is.numeric(rate) || length(rate) > 3){
    stop("Incorrect value for [.01, .99] bounded arguments", call. = FALSE)
  }

  if(!is.logical(two.tailed) || length(two.tailed) > 1){
    stop("Non-logical value for 'two.tailed'", call. = FALSE)
  }

  if(any(n3 - g3 < 3)){
    stop("Insufficient sample size, increase 'n3'", call. = FALSE)
  }


  df <- n3 - g3 - 2 - (n3 - 2) * (1 - n4)
  d <- .d(p, k1, k2, dists, rhots)
  sse <- sqrt(d * (rho3 * (1 - r23) / (p * (1 - p) * n4 * n3) +
                rho2 * (1 - r22) / (p * (1 - p) * n4 * n3 * n2) +
               (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n4 * n3 * n2 * n1)))

  mdes <- .mdes(power, alpha, sse, df, two.tailed)
  colnames(mdes) <- c("mdes", paste0(100 * (1 - round(alpha, 2)), "%lcl"),
                      paste0(100 * (1 - round(alpha, 2)), "%ucl"))
  mdes.out <- list(parms = list(power = power, alpha = alpha, two.tailed = two.tailed,
                                rhots = rhots, k1 = k1, k2 = k2, dists = dists,
                                rho2 = rho2, rho3 = rho3, r21 = r21, r22 = r22, r23 = r23,
                                g3 = g3, p = p, n1 = n1, n2 = n2, n3 = n3, n4 = n4),
                   df = df,
                   sse = sse,
                   mdes = mdes)
  class(mdes.out) <- c("mdes", "bcrd4f3")
  .summary.mdes(mdes.out)
  return(invisible(mdes.out))
}

# example
# mdes.crd3r3(rho3 = .06, rho2 = .17, n1 = 15, n2 = 3, n3 = 60)
# mdes.bcrd4f3(rho3 = .06, rho2 = .17, n1 = 15, n2 = 3, n3 = 60, n4 = 1)


power.bcrd4f3 <- function(es = .25, alpha = .05, two.tailed = TRUE,
                         rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                         rho2, rho3, r21 = 0, r22 = 0, r23 = 0, g3 = 0,
                         p = .50, n1, n2, n3, n4) {

  ss <- c(n1, n2, n3, n4, g3)
  if(any(ss < 0) || !is.numeric(ss) || length(ss) > 5){
    stop("Incorrect value for sample size or number of covariates", call. = FALSE)
  }

  var <- c(r21, r22, r23, rho2, rho3)
  if(any(var < 0) || any(var > 1) || !is.numeric(var) || length(var) > 5){
    stop("Incorrect value for [0, 1] bounded arguments", call. = FALSE)
  }

  rate <- c(p, alpha)
  if(any(rate < .01) || any(rate > .99) || !is.numeric(rate) || length(rate) > 2){
    stop("Incorrect value for [.01, .99] bounded arguments", call. = FALSE)
  }

  if(!is.logical(two.tailed) || length(two.tailed) > 1){
    stop("Non-logical value for 'two.tailed'", call. = FALSE)
  }

  if(any(n3 - g3 < 3)){
    stop("Insufficient sample size, increase 'n3'", call. = FALSE)
  }

  if(any(es <= 0) || !is.numeric(es) || length(es) > 1){
    stop("Incorrect value for 'es'", call. = FALSE)
  }

  df <- n3 - g3 - 2 - (n3 - 2) * (1 - n4)
  d <- .d(p, k1, k2, dists, rhots)
  sse <- sqrt(d * (rho3 * (1 - r23) / (p * (1 - p) * n4 * n3) +
                     rho2 * (1 - r22) / (p * (1 - p) * n4 * n3 * n2) +
                     (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n4 * n3 * n2 * n1)))

  power <- .power(es, alpha, sse, df, two.tailed)
  power.out <-  list(parms = list(es = es, alpha = alpha, two.tailed = two.tailed,
                                  rhots = rhots, k1 = k1, k2 = k2, dists = dists,
                                  rho2 = rho2, rho3 = rho3, r21 = r21, r22 = r22, r23 = r23,
                                  g3 = g3, p = p, n1 = n1, n2 = n2, n3 = n3, n4 = n4),
                     df = df,
                     sse = sse,
                     power = power)
  class(power.out) <- c("power", "bcrd4f3")
  .summary.power(power.out)
  return(invisible(power.out))
}
# example
# power.crd3r3(es = .446, rho3 = .06, rho2 = .17, n1 = 15, n2 = 3, n3 = 60)
# power.bcrd4f3(es = .446, rho3 = .06, rho2 = .17, n1 = 15, n2 = 3, n3 = 60, n4 = 1)

cosa.bcrd4f3 <- function(cn1 = 0, cn2 = 0, cn3 = 0, cn4 = 0, cost = NULL,
                        n1 = NULL, n2 = NULL, n3 = NULL, n4 = NULL,
                        p = NULL, n0 = c(10, 3, 100 + g3, 5), p0 = .499,
                        constrain = "power", round = TRUE, max.power = FALSE,
                        local.solver = c("LBFGS", "SLSQP", "MMA", "COBYLA"),
                        rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                        power = .80, es = .25, alpha = .05, two.tailed = TRUE,
                        rho2, rho3, g3 = 0, r21 = 0, r22 = 0, r23 = 0) {

  ss <- c(n1, n2, n3, n4, g3)
  if(any(ss < 0) || !is.numeric(ss) || length(ss) > 9){
    stop("Incorrect value for sample size or number of covariates", call. = FALSE)
  }

  var <- c(r21, r22, r23, rho2, rho3)
  if(any(var < 0) || any(var > 1) || !is.numeric(var) || length(var) > 5){
    stop("Incorrect value for [0, 1] bounded arguments", call. = FALSE)
  }

  rate <- c(alpha, power)
  if(any(rate < .01) || any(rate > .99) || !is.numeric(rate) || length(rate) > 2){
    stop("Incorrect value for [.01, .99] bounded arguments", call. = FALSE)
  }

  if(!is.logical(two.tailed) || length(two.tailed) > 1){
    stop("Non-logical value for 'two.tailed'", call. = FALSE)
  }

  if(!is.logical(max.power) || length(max.power) > 1){
    stop("Non-logical value for 'max.power'", call. = FALSE)
  }

  if(any(n3 - g3 < 3)){
    stop("Insufficient sample size, increase 'n3'", call. = FALSE)
  }

  if(any(es <= 0) || !is.numeric(es) || length(es) > 1){
    stop("Incorrect value for 'es'", call. = FALSE)
  }

  fun <- "cosa.bcrd4f3"
  lb <- c(1, 1, g3 + 3, 1)

  .df <- quote(n4 * (n3 - 2) - g3)
  .sse <- quote(sqrt(d * (rho3 * (1 - r23) / (p * (1 - p) * n4 * n3) +
                            rho2 * (1 - r22) / (p * (1 - p) * n4 * n3 * n2) +
                            (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n4 * n3 * n2 * n1))))
  .cost <- quote(n4 * cn4 +
                   n4 * n3 * (cn3[2] + p * (cn3[1] - cn3[2])) +
                   n4 * n3 * n2 * (cn2[2] + p * (cn2[1] - cn2[2])) +
                   n4 * n3 * n2 * n1 * (cn1[2] + p * (cn1[1] - cn1[2])))

  # NOTE: numerical derivatives for equality constraint (on power) are not precise
  # possibly due to numerical instability resulting from small step size (machine precision)
  # or round-off errors specific to designs with fixed effects
  .var.jacob <- expression(
    c(
      -d * (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n2 * n3 * n4 * n1^2),

      -d * rho2 * (1 - r22) / (p * (1 - p) * n2^2 * n3 * n4) -
        d * (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n2^2 * n3 * n4 * n1),

      -d * rho3 * (1 - r23) / (p * (1 - p) * n3^2 * n4) -
        d * rho2 * (1 - r22) / (p * (1 - p) * n2 * n3^2 * n4) -
        d * (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n2 * n3^2 * n4 * n1),

      -d * rho3 * (1 - r23) / (p * (1 - p) * n3 * n4^2) -
        d * rho2 * (1 - r22) / (p * (1 - p) * n2 * n3 * n4^2) -
        d * (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n2 * n3 * n4^2 * n1),

      -(1 - 2 * p) * d * rho3 * (1 - r23) / ((1 - p)^2 * p^2 * n3 * n4) -
        (1 - 2 * p) * d * rho2 * (1 - r22) / ((1 - p)^2 * p^2 * n2 * n3 * n4) -
        (1 - 2 * p) * d * (1 - rho3 - rho2) * (1 - r21) / ((1 - p)^2 * p^2 * n2 * n3 * n4 * n1)
    )
  )

  .cost.jacob <- expression(
    c(
      n4 * n3 * n2 * (p * cn1[1] + (1 - p) * cn1[2]),

      n4 * n3 * (p * cn2[1] + (1 - p) * cn2[2]) +
        n4 * n3 * n1 * (p * cn1[1] + (1 - p) * cn1[2]),

      n4 * (p * cn3[1] + (1 - p) * cn3[2]) +
        n4 * n2 * (p * cn2[1] + (1 - p) * cn2[2]) +
        n4 * n2 * n1 * (p * cn1[1] + (1 - p) * cn1[2]),

      cn4 +
        n3 * (p * cn3[1] + (1 - p) * cn3[2]) +
        n3 * n2 * (p * cn2[1] + (1 - p) * cn2[2]) +
        n3 * n2 * n1 * (p * cn1[1] + (1 - p) * cn1[2]),

      n4 * n3 * (cn3[1] - cn3[2]) +
        n4 * n3 * n2 * (cn2[1] - cn2[2]) +
        n4 * n3 * n2 * n1 * (cn1[1] - cn1[2])
    )
  )

  cosa <- .cosa(cn1 = cn1, cn2 = cn2, cn3 = cn3, cn4 = cn4, cost = cost,
                constrain = constrain, round = round,
                max.power = max.power, local.solver = local.solver,
                power = power, es = es, alpha = alpha, two.tailed = two.tailed,
                rhots = rhots, k1 = k1, k2 = k2, dists = dists,
                rho2 = rho2, rho3 = rho3, r21 = r21, r22 = r22, r23 = r23,
                g3 = g3, p0 = p0, p = p, n0 = n0, n1 = n1, n2 = n2, n3 = n3, n4 = n4)
  cosa.out <- list(parms = list(cn1 = cn1, cn2 = cn2, cn3 = cn3, cn4 = cn4, cost = cost,
                                constrain = constrain, round = round,
                                max.power = max.power, local.solver = local.solver,
                                power = power, es = es, alpha = alpha, two.tailed = two.tailed,
                                rhots = rhots, k1 = k1, k2 = k2, dists = dists,
                                rho2 = rho2, rho3 = rho3, r21 = r21, r22 = r22, r23 = r23,
                                g3 = g3, p0 = p0, p = p, n0 = n0, n1 = n1, n2 = n2, n3 = n3, n4 = n4),
                   cosa = cosa)
  class(cosa.out) <- c("cosa", "bcrd4f3")
  .summary.cosa(cosa.out)
  return(invisible(cosa.out))
}
# examples
# unconstrained
# cosa.crd3r3(constrain = "power", rho2 = .20, rho3 = .10, p = .35)
# cosa.bcrd4f3(constrain = "power", rho2 = .20, rho3 = .10, p = .35, n4 = 1)
# constrained
# cosa.crd3r3(cn1 = c(10,5), cn2 = c(20,10), cn3 = c(50,25), cost = 10000,
#             constrain = "cost", rho2 = .20, rho3 = .10, p = .35)
# cosa.bcrd4f3(cn1 = c(10,5), cn2 = c(20,10), cn3 = c(50,25), cost = 10000,
#            constrain = "cost", rho2 = .20, rho3 = .10, p = .35, n4 = 1)
