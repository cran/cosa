mdes.bird3r1 <- function(score = NULL, order = 2, rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                         power = .80, alpha = .05, two.tailed = TRUE, df = n3 - g3 - 1,
                         rho2, rho3, omega2, omega3, r21 = 0, r2t2 = 0, r2t3 = 0, g3 = 0,
                         rate.tp = 1, rate.cc = 0, p = .50, n1, n2, n3) {

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  if(df < 1) stop("Insufficient degrees of freedom", call. = FALSE)
  if(!is.null(score) && order == 0) warning("Ignoring information from score object \n", call. = FALSE)

  if(!is.null(rhots)) {
    if(rhots == 0) {
      if(order != 0) {
        order <- 0
        warning("'order' argument is ignored, 'order = 0' because 'rhots = 0'", call. = FALSE)
      }
      warning("'rhots' argument will be removed in the future,
               for corresponding random assignment designs use 'order = 0' instead", call. = FALSE)
    } else {
      stop("'rhots' argument will be removed in the future, arbitrary correlations are not allowed,
              use inspect.score() function instead", call. = FALSE)
    }
  }

  if(is.null(score)) {
    capture.output({
      score <- inspect.score(p = p, k1 = k1, k2 = k2, dists = dists)
    })
  } else {
    if("p" %in% names(user.parms)) {
      warning("Using 'p' from 'score' object, ignoring 'p' in the function call", call. = FALSE)
    }
    if(!inherits(score, "score")) {
      stop("'score' object should inherit 'score' class", call. = FALSE)
    }
  }

  p <- score$p
  ifelse(order == 2,
         d <- score$d2,
         ifelse(order == 1,
                d <- score$d1,
                d <- 1))

  if(order == 0) {
    score$parms$dists <- "uniform"
    score$parms$k1 <- 0
    score$parms$k2 <- 1
  }

  sse <- (1/(rate.tp - rate.cc)) * sqrt(rho3 * omega3 * (1 - r2t3) / n3 +
                rho2 * omega2 * (1 - r2t2) / (n3 * n2) +
                d * (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n3 * n2 * n1))

  mdes <- .mdes(power, alpha, sse, df, two.tailed)
  colnames(mdes) <- c("mdes", paste0(100 * (1 - round(alpha, 2)), "%lcl"),
                      paste0(100 * (1 - round(alpha, 2)), "%ucl"))
  mdes.out <- list(parms = list(order = order, dists = score$parms$dists,
                                k1 = score$parms$k1, k2 = score$parms$k2,
                                order = order,
                                power = power, alpha = alpha, two.tailed = two.tailed,
                                rho2 = rho2, rho3 = rho3, omega2 = omega2, omega3 = omega3,
                                r21 = r21, r2t2 = r2t2, r2t3 = r2t3,
                                g3 = g3, rate.tp = rate.tp, rate.cc = rate.cc,
                                p = p, n1 = n1, n2 = n2, n3 = n3),
                   df = df,
                   sse = sse,
                   mdes = mdes)
  class(mdes.out) <- c("mdes", "bird3r1")
  .summary.mdes(mdes.out)
  return(invisible(mdes.out))
}
mdes.bird3 <- mdes.bird3r1

power.bird3r1 <- function(score = NULL, order = 2, rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                          es = .25, alpha = .05, two.tailed = TRUE, df = n3 - g3 - 1,
                          rho2, rho3, omega2, omega3, r21 = 0, r2t2 = 0, r2t3 = 0, g3 = 0,
                          rate.tp = 1, rate.cc = 0, p = .50, n1, n2, n3) {

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  if(df < 1) stop("Insufficient degrees of freedom", call. = FALSE)
  if(!is.null(score) && order == 0) warning("Ignoring information from score object \n", call. = FALSE)

  if(!is.null(rhots)) {
    if(rhots == 0) {
      if(order != 0) {
        order <- 0
        warning("'order' argument is ignored, 'order = 0' because 'rhots = 0'", call. = FALSE)
      }
      warning("'rhots' argument will be removed in the future,
               for corresponding random assignment designs use 'order = 0' instead", call. = FALSE)
    } else {
      stop("'rhots' argument will be removed in the future, arbitrary correlations are not allowed,
              use inspect.score() function instead", call. = FALSE)
    }
  }

  if(is.null(score)) {
    capture.output({
      score <- inspect.score(p = p, k1 = k1, k2 = k2, dists = dists)
    })
  } else {
    if("p" %in% names(user.parms)) {
      warning("Using 'p' from 'score' object, ignoring 'p' in the function call", call. = FALSE)
    }
    if(!inherits(score, "score")) {
      stop("'score' object should inherit 'score' class", call. = FALSE)
    }
  }

  p <- score$p
  ifelse(order == 2,
         d <- score$d2,
         ifelse(order == 1,
                d <- score$d1,
                d <- 1))

  if(order == 0) {
    score$parms$dists <- "uniform"
    score$parms$k1 <- 0
    score$parms$k2 <- 1
  }

  sse <- (1/(rate.tp - rate.cc)) * sqrt(rho3 * omega3 * (1 - r2t3) / n3 +
                rho2 * omega2 * (1 - r2t2) / (n3 * n2) +
                d * (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n3 * n2 * n1))

  power <- .power(es, alpha, sse, df, two.tailed)
  power.out <-  list(parms = list(order = order, dists = score$parms$dists,
                                  k1 = score$parms$k1, k2 = score$parms$k2,
                                  es = es, alpha = alpha, two.tailed = two.tailed,
                                  rho2 = rho2, rho3 = rho3, omega2 = omega2, omega3 = omega3,
                                  r21 = r21, r2t2 = r2t2, r2t3 = r2t3,
                                  g3 = g3, rate.tp = rate.tp, rate.cc = rate.cc,
                                  p = p, n1 = n1, n2 = n2, n3 = n3),
                     df = df,
                     sse = sse,
                     power = power)
  names(power) <- "power"
  class(power.out) <- c("power", "bird3r1")
  .summary.power(power.out)
  return(invisible(power.out))
}
power.bird3 <-power.bird3r1

cosa.bird3r1 <- function(score = NULL, order = 2, rhots = NULL,
                         k1 = -6, k2 = 6, dists = "normal",
                         cn1 = 0, cn2 = 0, cn3 = 0, cost = NULL,
                         n1 = NULL, n2 = NULL, n3 = NULL, p = NULL,
                         n0 = c(10, 3, 100 + g3), p0 = .499,
                         constrain = "power", round = TRUE, max.power = FALSE,
                         local.solver = c("LBFGS", "SLSQP"),
                         power = .80, es = .25, alpha = .05, two.tailed = TRUE,
                         rho2, rho3, omega2, omega3,
                         g3 = 0, r21 = 0, r2t2 = 0, r2t3 = 0) {

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  if(!is.null(rhots)) {
    if(rhots == 0) {
      if(order != 0) {
        order <- 0
        warning("'order' argument is ignored, 'order = 0' because 'rhots = 0'", call. = FALSE)
      }
      warning("'rhots' argument will be removed in the future,
               for corresponding random assignment designs use 'order = 0' instead", call. = FALSE)
    } else {
      stop("'rhots' argument will be removed in the future, arbitrary correlations are not allowed,
              use inspect.score() function instead", call. = FALSE)
    }
  }

  if(!is.null(p)) {
    if(is.null(score)) {
      capture.output({
        score <- inspect.score(p = p, k1 = k1, k2 = k2, dists = dists)
      })
    } else {
      warning("Using 'p' from 'score' object, ignoring 'p' in the function call", call. = FALSE)
    }
    p <- score$p
    ifelse(order == 2,
           d <- score$d2,
           ifelse(order == 1,
                  d <- score$d1,
                  d <- 1))
  } else {
    if(order > 0) {
      if(is.null(score)) {
        stop("'order > 0' & 'p = NULL', p' cannot be NULL in regression discontinuity designs", call. = FALSE)
      } else {
        p <- score$p
        ifelse(order == 2,
               d <- score$d2,
               d <- score$d1)
      }
    } else if(order == 0){
      if(!is.null(score)) warning("Ignoring information from score object \n", call. = FALSE)
      d <- 1
      score$parms$dists <- "uniform"
      score$parms$k1 <- 0
      score$parms$k2 <- 1
    }
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

  .var.jacob <- expression(
    c(
      -d * (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n2 * n3 * n1^2),

      -rho2 * omega2 * (1 - r2t2) / (n2^2 * n3) -
        d * (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n2^2 * n3 * n1),

      -rho3 * omega3 * (1 - r2t3) / n3^2 -
        rho2 * omega2 * (1 - r2t2) / (n2 * n3^2) -
        d * (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n2 * n3^2 * n1),

      -(1 - 2 * p) * d * (1 - rho3 - rho2) * (1 - r21) / ((1 - p)^2 * p^2 * n2 * n3 * n1)
    )
  )

  .cost.jacob <- expression(
    c(
      n3 * n2 * (p * cn1[1] + (1 - p) * cn1[2]),

      n3 * cn2 +
        n3 * n1 * (p * cn1[1] + (1 - p) * cn1[2]),

      cn3 +
        n2 * cn2 +
        n2 * n1 * (p * cn1[1] + (1 - p) * cn1[2]),

      n3 * n2 * n1 * (cn1[1] - cn1[2])
    )
  )


  cosa <- .cosa(cn1 = cn1, cn2 = cn2, cn3 = cn3, cost = cost,
                constrain = constrain, round = round,
                max.power = max.power, local.solver = local.solver,
                power = power, es = es, alpha = alpha, two.tailed = two.tailed,
                rho2 = rho2, rho3 = rho3,
                omega2 = omega2, omega3 = omega3,
                r21 = r21, r2t2 = r2t2, r2t3 = r2t3,
                g3 = g3, p0 = p0, p = p, n0 = n0,
                n1 = n1, n2 = n2, n3 = n3)
  cosa.out <- list(parms = list(cn1 = cn1, cn2 = cn2, cn3 = cn3, cost = cost,
                                constrain = constrain, round = round,
                                max.power = max.power, local.solver = local.solver,
                                power = power, es = es, alpha = alpha, two.tailed = two.tailed,
                                order = order, dists = score$parms$dists,
                                k1 = score$parms$k1, k2 = score$parms$k2,
                                rho2 = rho2, rho3 = rho3,
                                omega2 = omega2, omega3 = omega3,
                                r21 = r21, r2t2 = r2t2, r2t3 = r2t3,
                                g3 = g3, p0 = p0, p = p, n0 = n0,
                                n1 = n1, n2 = n2, n3 = n3),
                   cosa = cosa)
  class(cosa.out) <- c("cosa", "bird3r1")
  .summary.cosa(cosa.out)
  return(invisible(cosa.out))
}
cosa.bird3 <- cosa.bird3r1
