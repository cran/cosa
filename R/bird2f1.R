mdes.bird2f1 <- function(power = .80, alpha = .05, two.tailed = TRUE,
                        rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                        r21 = 0, g1 = 0, p = .50, n1, n2 = 1) {

  ss <- c(n1, n2, g1)
  if(any(ss < 0) || !is.numeric(ss) || length(ss) > 3){
    stop("Incorrect value for sample size or number of covariates", call. = FALSE)
  }

  if(any(r21 < 0) || any(r21 > 1) || !is.numeric(r21) || length(r21) > 1){
    stop("Incorrect value for [0, 1] bounded arguments", call. = FALSE)
  }

  rate <- c(p, power, alpha)
  if(any(rate < .01) || any(rate > .99) || !is.numeric(rate) || length(rate) > 3){
    stop("Incorrect value for [.01, .99] bounded arguments", call. = FALSE)
  }

  if(!is.logical(two.tailed) || length(two.tailed) > 1){
    stop("Non-logical value for 'two.tailed'", call. = FALSE)
  }

  if(any(n1 - g1 < 3)){
    stop("Insufficient sample size, increase 'n1'", call. = FALSE)
  }

  user.call <- as.list(sys.call())
  names.user.call <- names(user.call)
  if("r2" %in%  names.user.call) {
    r21 <- user.call$r2
    warning("'r2' is renamed as 'r21' and will be removed from the next version", call. = FALSE)
  }
  if("g" %in%  names.user.call) {
    g1 <- user.call$g
    warning("'g' is renamed as 'g1' and will be removed from the next version", call. = FALSE)
  }
  if("n" %in%  names.user.call) {
    n1 <- user.call$n
    warning("'n' is renamed as 'n1' and will be removed from the next version", call. = FALSE)
  }

  df <- n1 - g1 - 2 - (n1 - 2) * (1 - n2)
  d <- .d(p, k1, k2, dists, rhots)
  sse <- sqrt(d * (1 - r21) / (p * (1 - p) * n1 * n2))

  mdes <- .mdes(power, alpha, sse, df, two.tailed)
  colnames(mdes) <- c("mdes", paste0(100 * (1 - round(alpha, 2)), "%lcl"),
                      paste0(100 * (1 - round(alpha, 2)), "%ucl"))
  mdes.out <- list(parms = list(power = power, alpha = alpha, two.tailed = two.tailed,
                                rhots = rhots, k1 = k1, k2 = k2, dists = dists,
                                r21 = r21, g1 = g1, p = p, n1 = n1, n2 = n2),
                   df = df,
                   sse = sse,
                   mdes = mdes)
  print(round(mdes, 3))
  class(mdes.out) <- c("mdes", "bird2f1")
  return(invisible(mdes.out))
}

# example
# mdes.bird2f1(n1 = 400, g1 = 3, r21 = .50)

power.bird2f1 <- function(es = .25, alpha = .05, two.tailed = TRUE,
                         rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                         r21 = 0, g1 = 0,
                         p = .50, n1, n2 = 1) {

  ss <- c(n1, n2, g1)
  if(any(ss < 0) || !is.numeric(ss) || length(ss) > 3){
    stop("Incorrect value for sample size or number of covariates", call. = FALSE)
  }

  if(any(r21 < 0) || any(r21 > 1) || !is.numeric(r21) || length(r21) > 1){
    stop("Incorrect value for [0, 1] bounded arguments", call. = FALSE)
  }

  rate <- c(p, alpha)
  if(any(rate < .01) || any(rate > .99) || !is.numeric(rate) || length(rate) > 2){
    stop("Incorrect value for [.01, .99] bounded arguments", call. = FALSE)
  }

  if(!is.logical(two.tailed) || length(two.tailed) > 1){
    stop("Non-logical value for 'two.tailed'", call. = FALSE)
  }

  if(any(n1 - g1 < 3)){
    stop("Insufficient sample size, increase 'n1'", call. = FALSE)
  }

  if(any(es <= 0) || !is.numeric(es) || length(es) > 1){
    stop("Incorrect value for 'es'", call. = FALSE)
  }

  user.call <- as.list(sys.call())
  names.user.call <- names(user.call)
  if("r2" %in%  names.user.call) {
    r21 <- user.call$r2
    warning("'r2' is renamed as 'r21' and will be removed from the next version", call. = FALSE)
  }
  if("g" %in%  names.user.call) {
    g1 <- user.call$g
    warning("'g' is renamed as 'g1' and will be removed from the next version", call. = FALSE)
  }
  if("n" %in%  names.user.call) {
    n1 <- user.call$n
    warning("'n' is renamed as 'n1' and will be removed from the next version", call. = FALSE)
  }

  df <- n1 - g1 - 2 - (n1 - 2) * (1 - n2)
  d <- .d(p, k1, k2, dists, rhots)
  sse <- sqrt(d * (1 - r21) / (p * (1 - p) * n1 * n2))

  power <- .power(es, alpha, sse, df, two.tailed)
  power.out <-  list(parms = list(es = es, alpha = alpha, two.tailed = two.tailed,
                                  rhots = rhots, k1 = k1, k2 = k2, dists = dists,
                                  r21 = r21, g1 = g1, p = p, n1 = n1, n2 = n2),
                     df = df,
                     sse = sse,
                     power = power)
  names(power) <- "power"
  print(round(power, 3))
  class(power.out) <- c("power", "bird2f1")
  return(invisible(power.out))
}
# example
# power.bird2f1(es = 0.466, n1 = 400, g1 = 3)

cosa.bird2f1 <- function(cn1 = 0, cn2 = 0, cost = NULL,
                        n1 = NULL, n2 = NULL, p = NULL, n0 = c(400 + g1, 5), p0 = .50,
                        constrain = "power", round = TRUE,
                        local.solver = c("LBFGS", "SLSQP", "MMA", "COBYLA"),
                        rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                        power = .80, es = .25, alpha = .05, two.tailed = TRUE,
                        g1 = 0, r21 = 0) {

  ss <- c(n1, n2, g1)
  if(any(ss < 0) || !is.numeric(ss) || length(ss) > 5){
    stop("Incorrect value for sample size or number of covariates", call. = FALSE)
  }

  if(any(r21 < 0) || any(r21 > 1) || !is.numeric(r21) || length(r21) > 1){
    stop("Incorrect value for [0, 1] bounded arguments", call. = FALSE)
  }

  rate <- c(alpha, power)
  if(any(rate < .01) || any(rate > .99) || !is.numeric(rate) || length(rate) > 2){
    stop("Incorrect value for [.01, .99] bounded arguments", call. = FALSE)
  }

  if(!is.logical(two.tailed) || length(two.tailed) > 1){
    stop("Non-logical value for 'two.tailed'", call. = FALSE)
  }

  if(any(n1 - g1 < 3)){
    stop("Insufficient sample size, increase 'n1'", call. = FALSE)
  }

  if(any(es <= 0) || !is.numeric(es) || length(es) > 1){
    stop("Incorrect value for 'es'", call. = FALSE)
  }

  fun <- "cosa.bird2f1"
  lb <- c(g1 + 3, 1)

  .df <- quote(n1 - g1 - 2 - (n1 - 2) * (1 - n2))
  .sse <- quote(sqrt(d * (1 - r21) / (p * (1 - p) * n2 * n1)))
  .cost <- quote(n2 * cn2 +
                   n2 * n1 * (cn1[2] + p * (cn1[1] - cn1[2])))

  cosa <- .cosa(cn1 = cn1, cn2 = cn2, cost = cost,
                constrain = constrain, round = round, local.solver = local.solver,
                power = power, es = es, alpha = alpha, two.tailed = two.tailed,
                rhots = rhots, k1 = k1, k2 = k2, dists = dists,
                r21 = r21, g1 = g1, p0 = p0, p = p, n0 = n0, n1 = n1, n2 = n2)
  cosa.out <- list(parms = list(cn1 = cn1, cn2 = cn2, cost = cost,
                                constrain = constrain, round = round, local.solver = local.solver,
                                power = power, es = es, alpha = alpha, two.tailed = two.tailed,
                                rhots = rhots, k1 = k1, k2 = k2, dists = dists,
                                r21 = r21, g1 = g1, p0 = p0, p = p, n0 = n0, n1 = n1, n2 = n2),
                   cosa = cosa)
  print(round(cosa, 3))
  class(cosa.out) <- c("cosa", "bird2f1")
  return(invisible(cosa.out))
}
# examples
# unconstrained
# cosa.ird1r1(rhots = 0)
# cosa.bird2f1(rhots = 0, n2 = 1)
# constrained
# cosa.ird1r1(cn1 = c(10,5), cost = 1000, constrain = "cost", p = .50)
# cosa.bird2f1(cn1 = c(10,5), cn2 = 0, n2 = 1, cost = 1000, constrain = "cost", p = .50)
# cosa.ird1r1(rhots = 0, cn1 = c(10,5), cost = 1000, constrain = "cost")
# cosa.bird2f1(rhots = 0, cn1 = c(10,5), n2 = 1, cost = 1000, constrain = "cost")


