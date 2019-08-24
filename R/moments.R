tnorm.moment <- function(mu = 0, sigma = 1, k1 = -10, k2 = 10, order = 1, central = FALSE) {
  if(order < 0) stop("'order' should be 0 or greater", call. = FALSE)

  if(central) {
    k1 <- k1 - mu
    k2 <- k2 - mu
    mu <- 0
  }

  zk1 <- (k1 - mu) / sigma
  zk2 <- (k2 - mu) / sigma

  M0 <- 1
  M1 <- mu - sigma*(dnorm(zk2) - dnorm(zk1)) / (pnorm(zk2) - pnorm(zk1))

  if(order == 0) {
    return(1)
  } else if(order == 1) {
    return(M1)
  } else {
    M <- (order - 1) * sigma^2 * tnorm.moment(mu = mu, sigma = sigma, k1 = k1, k2 = k2, order = order - 2) +
      mu * tnorm.moment(mu = mu, sigma = sigma, k1 = k1, k2 = k2, order = order - 1) -
      (sigma*(k2^(order - 1) * dnorm(zk2) - k1^(order - 1) * dnorm(zk1)) / (pnorm(zk2) - pnorm(zk1)))
    return(M)
  }
}

# moments of uniform
unif.moment <- function(k1 = 0, k2 = 1, order = 1, central = FALSE) {
  if(order < 0) stop("'order' should be 0 or greater", call. = FALSE)
  if(central) {
    M1 <- (k1 + k2)/2
    k1 <- k1 - M1
    k2 <- k2 - M1
  }
  M <- (k2^(order + 1) - k1^(order + 1)) / ((order + 1) * (k2 - k1))
  return(M)
}

# empirical moments
emp.moment <- function(x, order = 1, central = FALSE, absolute = FALSE, na.rm = FALSE) {

  if(!is.numeric(x) | !is.vector(x)) stop("'x' should be a numeric vector", call. = FALSE)

  if(na.rm) {
    x <- x[!is.na(x)]
  }

  if(central) {
    x <- x - mean(x)
    if(absolute) {
      x <- abs(x)
    }
  }

  if(absolute) {
    x <- abs(x)
  }

  return(sum(x^order) / length(x))
}

# function to compute p, rhots, rhots2
inspect.score <- function(score = NULL, sim = FALSE, p = NULL, cutoff = NULL, treat.lower = FALSE,
                          mu = 0, sigma = 1, k1 = -1e10, k2 =  1e10,
                          dists = "normal", ndraw = 1000, nsim = 1000) {

  pars <- c(p, cutoff, mu, sigma, k1, k2, ndraw, nsim)
  if(any(!is.numeric(pars) || length(pars) > 8)){
    stop("Incorrect value for one of the numeric argument", call. = FALSE)
  }

  if(is.null(cutoff) & is.null(p)){
    stop("Specify one of the 'cutoff' or 'p' arguments", call. = FALSE)
  }

  if(!is.null(cutoff) & !is.null(p)){
    warning("Ignoring 'cutoff' using 'p'", call. = FALSE)
  }

  if(!is.null(cutoff) & is.null(p)) {
    ifelse(dists == "normal",
           q <- msm::ptnorm(cutoff, mean = mu, sd = sigma, lower = k1, upper = k2),
           q <- punif(cutoff, min = k1, max = k2))
    ifelse(treat.lower, p <- q, p <- 1 - q)
  }

  if(!is.null(p)) {
    ifelse(treat.lower, q <- p, q <- 1 - p)
    ifelse(dists == "normal",
           cutoff <- msm::qtnorm(q, mean = mu, sd = sigma, lower = k1, upper = k2),
           cutoff <- qunif(q, min = k1, max = k2))
  }

  if(is.null(score) & isTRUE(sim)) {
    temp <- matrix(NA, nrow = nsim, ncol = 8)
    for(i in 1:nsim) {
      ifelse(dists == "normal",
             score <- msm::rtnorm(ndraw, mean = mu, sd = sigma, lower = k1, upper = k2),
             score <- runif(ndraw, k1, k2))

      ifelse(isFALSE(treat.lower),
             treatment <- ifelse(score > cutoff, 1, 0),
             treatment <- ifelse(score < cutoff, 1, 0))
      temp[i,1] <- mean(treatment)
      temp[i,2] <- cor(treatment, score)
      temp[i,3] <- cor(treatment, score^2)
      temp[i,4] <- cor(score, score^2)
      temp[i,5] <- mean(score)
      temp[i,6] <- var(score)
      temp[i,7] <- min(score)
      temp[i,8] <- max(score)
    }
    p <- mean(temp[,1])
    rhots <- mean(temp[,2])
    rhots2 <- mean(temp[,3])
    rhoss2 <- mean(temp[,4])
    mu <- mean(temp[,5])
    sigma <- sqrt(sum(temp[,6]) / nsim)
    k1 <- mean(temp[,7])
    k2 <- mean(temp[,8])
  } else if(is.null(score) & isFALSE(sim)) {
    if(dists == "normal") {
      m1 <- tnorm.moment(mu = mu, sigma = sigma, k1 = k1, k2 = k2, central = FALSE, order = 1)
      m2 <- tnorm.moment(mu = mu, sigma = sigma, k1 = k1, k2 = k2, central = FALSE, order = 2)
      m3 <- tnorm.moment(mu = mu, sigma = sigma, k1 = k1, k2 = k2, central = FALSE, order = 3)
      m4 <- tnorm.moment(mu = mu, sigma = sigma, k1 = k1, k2 = k2, central = FALSE, order = 4)
      ifelse(treat.lower,
             m1z0 <- tnorm.moment(mu = mu, sigma = sigma, k1 = k1, k2 = cutoff, central = FALSE, order = 1),
             m1z0 <- tnorm.moment(mu = mu, sigma = sigma, k1 = cutoff, k2 = k2, central = FALSE, order = 1))
      ifelse(treat.lower,
             m2z0 <- tnorm.moment(mu = mu, sigma = sigma, k1 = k1, k2 = cutoff, central = FALSE, order = 2),
             m2z0 <- tnorm.moment(mu = mu, sigma = sigma, k1 = cutoff, k2 = k2, central = FALSE, order = 2))
      m4 <- tnorm.moment(mu = mu, sigma = sigma, k1 = k1, k2 = k2, central = FALSE, order = 4)
      rhots <- p*(m1z0 - m1) /  sqrt(p*(1-p)*(m2 - m1^2))
      rhots2 <- p*(m2z0 - m2) /  sqrt(p*(1-p)*(m4 - (m2)^2))
      rhoss2 <- (m3 - m1*m2) / sqrt((m2 - m1^2) * (m4 - m2^2))
      # mu <- m1 # biased
      # sigma <- sqrt(tnorm.moment(mu = mu, sigma = sigma, k1 = k1, k2 = k2, central = TRUE, order = 2)) # biased
      mu <- mu
      sigma <- sigma
    } else if(dists == "uniform") {
      m1 <- unif.moment(k1 = k1, k2 = k2, order = 1, central = FALSE)
      m2 <- unif.moment(k1 = k1, k2 = k2, order = 2, central = FALSE)
      m3 <- unif.moment(k1 = k1, k2 = k2, order = 3, central = FALSE)
      m4 <- unif.moment(k1 = k1, k2 = k2, order = 4, central = FALSE)
      ifelse(treat.lower,
             m1z0 <- unif.moment(k1 = k1, k2 = cutoff, order = 1, central = FALSE),
             m1z0 <- unif.moment(k1 = cutoff, k2 = k2, order = 1, central = FALSE))
      ifelse(treat.lower,
             m2z0 <- unif.moment(k1 = k1, k2 = cutoff, order = 2, central = FALSE),
             m2z0 <- unif.moment(k1 = cutoff, k2 = k2, order = 2, central = FALSE))
      m4 <- unif.moment(k1 = k1, k2 = k2, order = 4, central = FALSE)
      rhots <- p*(m1z0 - m1) /  sqrt(p*(1-p)*(m2 - m1^2))
      rhots2 <- p*(m2z0 - m2) /  sqrt(p*(1-p)*(m4 - (m2)^2))
      rhoss2 <- (m3 - m1*m2) / sqrt((m2 - m1^2) * (m4 - m2^2))
      mu <- m1
      sigma <- sqrt(unif.moment(k1 = k1, k2 = k2, order = 2, central = TRUE))
    } else if(dists == "empirical"){
      stop("Empirical score variable is missing", call. = FALSE)
    } else {
      stop("Distribution type not supported", call. = FALSE)
    }

  } else if(!is.null(score)){
    score <- score[score < k2 & score > k1]
    ifelse(isFALSE(treat.lower),
           treatment <- ifelse(score > cutoff, 1, 0),
           treatment <- ifelse(score < cutoff, 1, 0))
    cutoff <- quantile(score, q)
    p <- mean(treatment)
    rhots <- cor(treatment, score)
    rhots2 <- cor(treatment, score^2)
    rhoss2 <- cor(score, score^2)
    mu <- mean(score)
    sigma <- sd(score)
    k1 <- min(score)
    k2 <- max(score)
  }

  d1 <- 1 / (1 - rhots^2)
  d2 <- (1 - rhoss2^2) / (1 - rhots^2 - rhots2^2 - rhoss2^2 + 2 * rhots * rhots2 * rhoss2)

  if(!is.null(score)) dists <- "empirical"
  parms <- list(dists = dists, mu = mu, sigma = sigma, k1 = k1, k2 = k2)
  score.out <- list(parms = parms, cutoff = cutoff, treat.lower = treat.lower, p = p, rhots = rhots, rhots2 = rhots2, rhoss2 = rhoss2,
                    d1 = d1, d2 = d2)
  class(score.out) <- "score"
  cat("\nInspect score results \n------------------------------ \n")
  result <- cbind(p = p, cutoff = cutoff, rhots = rhots, rhots2 = rhots2, rhoss2 = rhoss2)
  print(round(result,2))
  return(invisible(score.out))
}
