mdes.ird1r1 <- function(power = .80, alpha = .05, two.tailed = TRUE,
                        rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                        r21 = 0, g1 = 0, p = .50, n1) {

  ss <- c(n1, g1)
  if(any(ss < 0) || !is.numeric(ss) || length(ss) > 2){
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

  df <- n1 - g1 - 2
  d <- .d(p, k1, k2, dists, rhots)
  sse <- sqrt(d * (1 - r21) / (p * (1 - p) * n1))

  mdes <- .mdes(power, alpha, sse, df, two.tailed)
  colnames(mdes) <- c("mdes", paste0(100 * (1 - round(alpha, 2)), "%lcl"),
                      paste0(100 * (1 - round(alpha, 2)), "%ucl"))
  mdes.out <- list(parms = list(power = power, alpha = alpha, two.tailed = two.tailed,
                                rhots = rhots, k1 = k1, k2 = k2, dists = dists,
                                r21 = r21, g1 = g1, p = p, n1 = n1),
                   df = df,
                   sse = sse,
                   mdes = mdes)
  class(mdes.out) <- c("mdes", "ird1r1")
  .summary.mdes(mdes.out)
  return(invisible(mdes.out))
}

# example
# mdes.ird1r1(n1 = 400, g1 = 3, r21 = .50)

power.ird1r1 <- function(es = .25, alpha = .05, two.tailed = TRUE,
                         rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
                         r21 = 0, g1 = 0,
                         p = .50, n1) {

  ss <- c(n1, g1)
  if(any(ss < 0) || !is.numeric(ss) || length(ss) > 2){
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

  df <- n1 - g1 - 2
  d <- .d(p, k1, k2, dists, rhots)
  sse <- sqrt(d * (1 - r21) / (p * (1 - p) * n1))

  power <- .power(es, alpha, sse, df, two.tailed)
  power.out <-  list(parms = list(es = es, alpha = alpha, two.tailed = two.tailed,
                                  rhots = rhots, k1 = k1, k2 = k2, dists = dists,
                                  r21 = r21, g1 = g1, p = p, n1 = n1),
                     df = df,
                     sse = sse,
                     power = power)
  class(power.out) <- c("power", "ird1r1")
  .summary.power(power.out)
  return(invisible(power.out))
}
# example
# power.ird1r1(es = 0.466, n1 = 400, g1 = 3)


