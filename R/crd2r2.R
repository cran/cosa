mdes.crd2r2 <- function(power = .80, alpha = .05, two.tailed = TRUE,
                        rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                        rho2, r21 = 0, r22 = 0, g2 = 0,
                        p = .50, n1, n2) {

  ss <- c(n1, n2, g2)
  if(any(ss < 0) || !is.numeric(ss) || length(ss) > 3){
    stop("Incorrect value for sample size or number of covariates", call. = FALSE)
  }

  var <- c(r21, r22, rho2)
  if(any(var < 0) || any(var > 1) || !is.numeric(var) || length(var) > 3){
    stop("Incorrect value for [0, 1] bounded arguments", call. = FALSE)
  }

  rate <- c(p, power, alpha)
  if(any(rate < .01) || any(rate > .99) || !is.numeric(rate) || length(rate) > 3){
    stop("Incorrect value for [.01, .99] bounded arguments", call. = FALSE)
  }

  if(!is.logical(two.tailed) || length(two.tailed) > 1){
    stop("Non-logical value for 'two.tailed'", call. = FALSE)
  }

  if(any(n2 - g2 < 3)){
    stop("Insufficient sample size, increase 'n2'", call. = FALSE)
  }

  df <- n2 - g2 - 2
  d <- .d(p, k1, k2, dists, rhots)
  sse <- sqrt(d * (rho2 * (1 - r22) / (p * (1 - p) * n2) +
               (1 - rho2) * (1 - r21) / (p * (1 - p) * n2 * n1)))

  mdes <- .mdes(power, alpha, sse, df, two.tailed)
  colnames(mdes) <- c("mdes", paste0(100 * (1 - round(alpha, 2)), "%lcl"),
                      paste0(100 * (1 - round(alpha, 2)), "%ucl"))
  mdes.out <- list(parms = list(power = power, alpha = alpha, two.tailed = two.tailed,
                                rhots = rhots, k1 = k1, k2 = k2, dists = dists,
                                rho2 = rho2, r21 = r21, r22 = r22,
                                g2 = g2, p = p, n1 = n1, n2 = n2),
                   df = df,
                   sse = sse,
                   mdes = mdes)
  class(mdes.out) <- c("mdes", "crd2r2")
  .summary.mdes(mdes.out)
  return(invisible(mdes.out))
}

# example
# mdes.crd2r2(rho2=.20, n1 = 4, n2 = 20)


power.crd2r2 <- function(es = .25, alpha = .05, two.tailed = TRUE,
                         rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                         rho2, r21 = 0, r22 = 0, g2 = 0,
                         p = .50, n1, n2) {

  ss <- c(n1, n2, g2)
  if(any(ss < 0) || !is.numeric(ss) || length(ss) > 3){
    stop("Incorrect value for sample size or number of covariates", call. = FALSE)
  }

  var <- c(r21, r22, rho2)
  if(any(var < 0) || any(var > 1) || !is.numeric(var) || length(var) > 3){
    stop("Incorrect value for [0, 1] bounded arguments", call. = FALSE)
  }

  rate <- c(p, alpha)
  if(any(rate < .01) || any(rate > .99) || !is.numeric(rate) || length(rate) > 2){
    stop("Incorrect value for [.01, .99] bounded arguments", call. = FALSE)
  }

  if(!is.logical(two.tailed) || length(two.tailed) > 1){
    stop("Non-logical value for 'two.tailed'", call. = FALSE)
  }

  if(any(n2 - g2 < 3)){
    stop("Insufficient sample size, increase 'n2'", call. = FALSE)
  }

  if(any(es <= 0) || !is.numeric(es) || length(es) > 1){
    stop("Incorrect value for 'es'", call. = FALSE)
  }

  df <- n2 - g2 - 2
  d <- .d(p, k1, k2, dists, rhots)
  sse <- sqrt(d * (rho2 * (1 - r22) / (p * (1 - p) * n2) +
                     (1 - rho2) * (1 - r21) / (p * (1 - p) * n2 * n1)))

  power <- .power(es, alpha, sse, df, two.tailed)
  power.out <-  list(parms = list(es = es, alpha = alpha, two.tailed = two.tailed,
                                  rhots = rhots, k1 = k1, k2 = k2, dists = dists,
                                  rho2 = rho2, r21 = r21, r22 = r22,
                                  g2 = g2, p = p, n1 = n1, n2 = n2),
                     df = df,
                     sse = sse,
                     power = power)
  class(power.out) <- c("power", "crd2r2")
  .summary.power(power.out)
  return(invisible(power.out))
}
# example
# power.crd2r2(es = 1.391, rho2 = .20, n1 = 4, n2 = 20)

cosa.crd2r2 <- function(cn1 = 0, cn2 = 0, cost = NULL,
                        n1 = NULL, n2 = NULL, p = NULL, n0 = c(10, 100 + g2), p0 = .499,
                        constrain = "power", round = TRUE, max.power = FALSE,
                        local.solver = c("LBFGS", "SLSQP"),
                        rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                        power = .80, es = .25, alpha = .05, two.tailed = TRUE,
                        rho2, g2 = 0, r21 = 0, r22 = 0) {

  ss <- c(n1, n2, g2)
  if(any(ss < 0) || !is.numeric(ss) || length(ss) > 5){
    stop("Incorrect value for sample size or number of covariates", call. = FALSE)
  }

  var <- c(r21, r22, rho2)
  if(any(var < 0) || any(var > 1) || !is.numeric(var) || length(var) > 3){
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

  if(any(n2 - g2 < 3)){
    stop("Insufficient sample size, increase 'n2'", call. = FALSE)
  }

  if(any(es <= 0) || !is.numeric(es) || length(es) > 1){
    stop("Incorrect value for 'es'", call. = FALSE)
  }

  fun <- "cosa.crd2r2"
  lb <- c(1, g2 + 3)

  .df <- quote(n2 - g2 - 2)
  .sse <- quote(sqrt(d * (rho2 * (1 - r22) / (p * (1 - p) * n2) +
                            (1 - rho2) * (1 - r21) / (p * (1 - p) * n2 * n1))))
  .cost <- quote(n2 * (cn2[2] + p * (cn2[1] - cn2[2])) +
                   n2 * n1 * (cn1[2] + p * (cn1[1] - cn1[2])))

  .var.jacob <- expression(
    c(
      -d * (1 - rho2) * (1 - r21) / (p * (1 - p) * n2 * n1^2),

      -d * rho2 * (1 - r22) / (p * (1 - p) * n2^2) -
        d * (1 - rho2) * (1 - r21) / (p * (1 - p) * n2^2 * n1),

      -(1 - 2 * p) * d * rho2 * (1 - r22) / ((1 - p)^2 * p^2 * n2) -
        (1 - 2 * p) * d * (1 - rho2) * (1 - r21) / ((1 - p)^2 * p^2 * n2 * n1)
    )
  )

  .cost.jacob <- expression(
    c(
      n2 * (p * cn1[1] + (1 - p) * cn1[2]),

      (p * cn2[1] + (1 - p) * cn2[2]) +
        n1 * (p * cn1[1] + (1 - p) * cn1[2]),

      n2 * (cn2[1] - cn2[2]) +
        n2 * n1 * (cn1[1] - cn1[2])
    )
  )

  cosa <- .cosa(cn1 = cn1, cn2 = cn2, cost = cost,
                constrain = constrain, round = round,
                max.power = max.power, local.solver = local.solver,
                power = power, es = es, alpha = alpha, two.tailed = two.tailed,
                rhots = rhots, k1 = k1, k2 = k2, dists = dists,
                rho2 = rho2, r21 = r21, r22 = r22,
                g2 = g2, p0 = p0, p = p, n0 = n0, n1 = n1, n2 = n2)
  cosa.out <- list(parms = list(cn1 = cn1, cn2 = cn2, cost = cost,
                                constrain = constrain, round = round,
                                max.power = max.power, local.solver = local.solver,
                                power = power, es = es, alpha = alpha, two.tailed = two.tailed,
                                rhots = rhots, k1 = k1, k2 = k2, dists = dists,
                                rho2 = rho2, r21 = r21, r22 = r22,
                                g2 = g2, p0 = p0, p = p, n0 = n0, n1 = n1, n2 = n2),
                   cosa = cosa)
  class(cosa.out) <- c("cosa", "crd2r2")
  .summary.cosa(cosa.out)
  return(invisible(cosa.out))
}
# examples
# unconstrained
# cosa.crd2r2(rhots = 0, constrain = "power", rho2 = .20, local.solver = "lbfgs")
# constrained
#
# cosa.crd2r2(cn1 = c(10,5), cn2 = c(20,10), cost = 5000, p = .35,
#            constrain = "cost", rho2 = .20)
