mdes.ird1r1 <- function(score = NULL, order = 2, rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                        power = .80, alpha = .05, two.tailed = TRUE, df = n1 - g1 - order - 2,
                        r21 = 0, g1 = 0, rate.tp = 1, rate.cc = 0, p = .50, n1) {

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

  sse <- (1/(rate.tp - rate.cc)) * sqrt(d * (1 - r21) / (p * (1 - p) * n1))

  mdes <- .mdes(power, alpha, sse, df, two.tailed)
  colnames(mdes) <- c("mdes", paste0(100 * (1 - round(alpha, 2)), "%lcl"),
                      paste0(100 * (1 - round(alpha, 2)), "%ucl"))
  mdes.out <- list(parms = list(order = order, power = power, alpha = alpha, two.tailed = two.tailed,
                                r21 = r21, g1 = g1, rate.tp = rate.tp, rate.cc = rate.cc, p = p, n1 = n1),
                   df = df,
                   sse = sse,
                   mdes = mdes)

  class(mdes.out) <- c("mdes", "ird1r1")
  .summary.mdes(mdes.out)
  return(invisible(mdes.out))
}
mdes.ird <- mdes.ird1r1

power.ird1r1 <- function(score = NULL, order = 2, rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                         es = .25, alpha = .05, two.tailed = TRUE, df = n1 - g1 - order - 2,
                         r21 = 0, g1 = 0, rate.tp = 1, rate.cc = 0, p = .50, n1) {


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

  sse <- (1/(rate.tp - rate.cc)) * sqrt(d * (1 - r21) / (p * (1 - p) * n1))

  power <- .power(es, alpha, sse, df, two.tailed)
  power.out <-  list(parms = list(order = order, es = es, alpha = alpha, two.tailed = two.tailed,
                                  r21 = r21, g1 = g1, rate.tp = rate.tp, rate.cc = rate.cc, p = p, n1 = n1),
                     df = df,
                     sse = sse,
                     power = power)
  class(power.out) <- c("power", "ird1r1")
  .summary.power(power.out)
  return(invisible(power.out))
}
power.ird <- power.ird1r1
