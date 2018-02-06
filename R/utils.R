# minimum detectable effect size
.mdes <- function(power, alpha, sse, df, two.tailed){
  t1 <- ifelse(two.tailed == TRUE, abs(qt(alpha / 2, df)), abs(qt(alpha, df)))
  t2 <- abs(qt(power, df))
  m <- ifelse(power >= 0.5, t1 + t2, t1 - t2)
  mdes <- m * sse
  lcl <- mdes * (1 - t1 / m)
  ucl <- mdes * (1 + t1 / m)
  mlu <- cbind(mdes, lcl, ucl)
  return(mlu)
}

# statistical power
.power <- function(es, alpha, sse, df, two.tailed){
  lambda <- es/sse
  power <- ifelse(two.tailed == FALSE,
                  1 - pt(qt(alpha, df, lower.tail = FALSE), df, lambda),
                  1 - pt(qt(alpha / 2, df, lower.tail = FALSE), df, lambda) +
                    pt(-qt(alpha / 2, df, lower.tail = FALSE), df, lambda))
  return(power)
}

# normal (Schochet, 2008, p.14)
# p: proportion of cases in treatment
.rhotsn <- function(p){
  rhots <- dnorm(qnorm(1 - p)) / sqrt(p * (1 - p))
  return(rhots)
}

# uniform (Schochet, 2008, p.14)
# p: proportion of cases in treatment
.rhotsu <- function(p){
  rhots <- sqrt(3 * p * (1 - p))
  return(rhots)
}

# truncated normal (Schochet, 2008, p.14)
# p: proportion of cases in treatment
# k1: left truncation point (in standard deviation units from full normal distribution mean)
# k2: right truncation point (in standard deviation units from full normal distribution mean)
.rhotstn <- function(k1, k2, p){
  c <- qnorm(p * pnorm(k1) + (1 - p) * pnorm(k2) )
  sigmas2 <- 1 - ((k2 * dnorm(k2) - k1 * dnorm(k1)) / (pnorm(k2) - pnorm(k1))) -
    ((dnorm(k2) - dnorm(k1)) / (pnorm(k2) - pnorm(k1)))^2
  rhots <- (p / sqrt(sigmas2 * p * (1 - p))) * ((dnorm(k2) - dnorm(k1)) / (pnorm(k2) - pnorm(k1)) -
                                                  (dnorm(k2) - dnorm(c)) / (pnorm(k2) - pnorm(c)))
  return(rhots)
}

# design effect
.d <- function(p, k1, k2, dists, rhots) {
  .de <- function(rhots) {
    d <- 1 / (1 - rhots^2)
    return(d)
  }
  if(is.null(rhots)) {
    if(dists == "normal") {
      d <- .de(.rhotstn(k1, k2, p))
    } else if(dists == "uniform") {
      if(k1 != -6 | k2 != 6) {
        warning("k1 and/or k2 will be ignored.", call. = FALSE)
      }
      d <- .de(.rhotsu(p))
    }
  } else if(!is.null(rhots)) {
    d <- .de(rhots)
  }
  return(d)
}

# constrained optimal sample allocation
.cosa <- function(cn1 = 0, cn2 = 0, cn3 = 0, cn4 = 0, cost = NULL,
                  n1 = NULL, n2 = NULL, n3 = NULL, n4 = NULL, p = NULL, n0, p0,
                  constrain, local.solver = c("LBFGS", "SLSQP", "MMA", "COBYLA"),
                  rhots, k1, k2, dists,
                  power, es, alpha, two.tailed,
                  rho2, rho3, rho4, omega2, omega3, omega4,
                  g1 = 0, g2 = 0, g3 = 0, g4 = 0,
                  r21 = 0, r22 = 0, r23 = 0, r24 = 0,
                  r2t2 = 0, r2t3 = 0, r2t4 = 0) {

  .df <- get(".df", parent.frame())
  .sse <- get(".sse", parent.frame())
  .cost <- get(".cost", parent.frame())
  lb <- get("lb", parent.frame())
  fun <- get("fun", parent.frame())

  d <- .d(p, k1, k2, dists, rhots)
  if(!is.null(rhots)) {
    if(rhots != 0 & is.null(p)) {
      stop("'p' cannot be NULL in regression discontinuity design. Specify 'rhots = 0' ", call. = FALSE)
    } else if(rhots == 0) {
      cat(" Results are equivalent to random assignment design \n")
    }
  } else if(is.null(rhots) & is.null(p)) {
    stop("'p' cannot be NULL in regression discontinuity design. Specify 'rhots = 0' ", call. = FALSE)
  }

  fun_parsed <- scan(text = fun, what = "character", sep=".", quiet = TRUE)
  rlevel <- as.numeric(substr(fun, nchar(fun), nchar(fun)))
  nlevels <- as.numeric(substr(fun, nchar(fun)-2, nchar(fun)-2))

  # equality constraints on cost
  .eq.cost <- function(ss){
    n1 <- ss[1]
    n2 <- ss[2]
    if(nlevels >= 3){
      n3 <- ss[3]
    }
    if(nlevels == 4){
      n4 <- ss[4]
    }
    p <- ss[nlevels + 1]
    return(eval(.cost) - cost)
  }

  # equality constraints on power
  .eq.power <- function(ss){
    n1 <- ss[1]
    n2 <- ss[2]
    if(nlevels >= 3){
      n3 <- ss[3]
    }
    if(nlevels == 4){
      n4 <- ss[4]
    }
    p <- ss[nlevels + 1]
    return(.power(es = es,  alpha = alpha,
                  sse = eval(.sse), df = eval(.df),
                  two.tailed = two.tailed) - power)
  }

  # equality constraints on mdes
  .eq.es <- function(ss){
    n1 <- ss[1]
    n2 <- ss[2]
    if(nlevels >= 3){
      n3 <- ss[3]
    }
    if(nlevels == 4){
      n4 <- ss[4]
    }
    p <- ss[nlevels + 1]
    return(.mdes(power = power, alpha = alpha,
                 sse = eval(.sse), df = eval(.df),
                 two.tailed = two.tailed)[1] - es)
  }

  #  mdes and power
  .mlu.pwr <- function(ss){
    n1 <- ss[1]
    n2 <- ss[2]
    if(nlevels >= 3){
      n3 <- ss[3]
    }
    if(nlevels == 4){
      n4 <- ss[4]
    }
    p <- ss[nlevels + 1]
    return(c(.mdes(power = power, alpha = alpha,
                   sse = eval(.sse), df = eval(.df),
                   two.tailed = two.tailed),
             .power(es = es, alpha = alpha,
                    sse = eval(.sse), df = eval(.df),
                    two.tailed = two.tailed)))
  }

  # minimize cost given power or mdes
  .min.cost <- function(ss){
    n1 <- ss[1]
    n2 <- ss[2]
    if(nlevels >= 3){
      n3 <- ss[3]
    }
    if(nlevels == 4){
      n4 <- ss[4]
    }
    p <- ss[nlevels + 1]
    return(eval(.cost))
  }

  # minimize treatment variance given cost
  .min.var <- function(ss){
    n1 <- ss[1]
    n2 <- ss[2]
    if(nlevels >= 3){
      n3 <- ss[3]
    }
    if(nlevels == 4){
      n4 <- ss[4]
    }
    p <- ss[nlevels + 1]
    return(eval(.sse)^2)
  }

  cost.list <- list(cn1, cn2, cn3, cn4)
  cost.names <- c("cn1", "cn2", "cn3", "cn4")
  cost.lengths <- unlist(lapply(cost.list, length))
  if(rlevel < nlevels & any(cost.lengths[(rlevel+1):nlevels] == 2)) {
    stop("Unequal cost applies to levels at or below randomization (or discontinuity) level", call.=FALSE)
  }

  if(any(cost.lengths[1:rlevel] > 2)) {
    stop("Marginal costs cannot have a length greater than two", call.=FALSE)
  }

  if(!is.null(cost)) {
    if(length(cost)>1) {
      stop("Total cost cannot have a length greater than one", call.=FALSE)
    }else if(cost<0) {
      stop(" Total cost cannot have a negative value", call.=FALSE)
    }
  }

  fn.constr <- switch(constrain,
                      "power" = .eq.power,
                      "es" = .eq.es,
                      "cost" = .eq.cost,
                      stop("Unknown constrain", call. = FALSE))
  fn.min <- switch(constrain,
                   "power" = .min.cost,
                   "es" = .min.cost,
                   "cost" = .min.var)

  if(rlevel >= 1 & length(cn1) == 1) {
    cn1 <- c(cn1, cn1)
  }
  if(rlevel >= 2 & length(cn2) == 1) {
    cn2 <- c(cn2, cn2)
  }
  if(rlevel >= 3 & length(cn3) == 1) {
    cn3 <- c(cn3,cn3)
  }
  if(rlevel == 4 & length(cn4) == 1) {
    cn4 <- c(cn4, cn4)
  }

  # constraints on one sample sizes and p
  p0  <- ifelse(!is.null(p), mean(p), p0)
  plb <- ifelse(!is.null(p), min(p), .15)
  pub <- ifelse(!is.null(p), max(p), .85)
  n10  <- ifelse(!is.null(n1), mean(n1), n0[1])
  n1lb <- ifelse(!is.null(n1), min(n1), lb[1])
  n1ub <- ifelse(!is.null(n1), max(n1), Inf)
  if(nlevels == 2){
    n20  <- ifelse(!is.null(n2), mean(n2), n0[2])
    n2lb <- ifelse(!is.null(n2), min(n2), lb[2])
    n2ub <- ifelse(!is.null(n2), max(n2), Inf)
    ss0  <- c(n10,n20,p0)
    sslb <- c(n1lb,n2lb,plb)
    ssub <- c(n1ub,n2ub,pub)
  }else if(nlevels == 3){
    n20  <- ifelse(!is.null(n2), mean(n2), n0[2])
    n2lb <- ifelse(!is.null(n2), min(n2), lb[2])
    n2ub <- ifelse(!is.null(n2), max(n2), Inf)
    n30  <- ifelse(!is.null(n3), mean(n3), n0[3])
    n3lb <- ifelse(!is.null(n3), min(n3), lb[3])
    n3ub <- ifelse(!is.null(n3), max(n3), Inf)
    ss0  <- c(n10,n20,n30,p0)
    sslb <- c(n1lb,n2lb,n3lb,plb)
    ssub <- c(n1ub,n2ub,n3ub,pub)
  }else if(nlevels == 4){
    n20  <- ifelse(!is.null(n2), mean(n2), n0[2])
    n2lb <- ifelse(!is.null(n2), min(n2), lb[2])
    n2ub <- ifelse(!is.null(n2), max(n2), Inf)
    n30  <- ifelse(!is.null(n3), mean(n3), n0[3])
    n3lb <- ifelse(!is.null(n3), min(n3), lb[3])
    n3ub <- ifelse(!is.null(n3), max(n3), Inf)
    n40  <- ifelse(!is.null(n4), mean(n4), n0[4])
    n4lb <- ifelse(!is.null(n4), min(n4), lb[4])
    n4ub <- ifelse(!is.null(n4), max(n4), Inf)
    ss0  <- c(n10,n20,n30,n40,p0)
    sslb <- c(n1lb,n2lb,n3lb,n4lb,plb)
    ssub <- c(n1ub,n2ub,n3ub,n4ub,pub)
  }

  # constrained optimal sample allocation
  #local.solver <- c("SLSQP", "LBFGS", "MMA", "COBYLA")
  i <- 1; conv <- FALSE
  while(i <= length(local.solver) & conv == FALSE){
    nlopt.ss <- nloptr::auglag(x0 = ss0, fn = fn.min, heq = fn.constr,
                               localsolver = local.solver[i], localtol = 1e-3,
                               lower = sslb, upper = ssub)
    if(nlopt.ss$par[nlevels] %% 1 != 0) {
      ss0[nlevels] <- sslb[nlevels] <- ssub[nlevels] <- round(nlopt.ss$par[nlevels])
      nlopt.ss <- nloptr::auglag(x0 = ss0, fn = fn.min, heq = fn.constr,
                                 localsolver = local.solver[i], localtol = 1e-3,
                                 lower = sslb, upper = ssub)
    }
    if(nlopt.ss$convergence < 0 | all(nlopt.ss$par == ss0) | any(nlopt.ss$par <= 0)) {
      conv <- FALSE
      cat(" Solution is not feasible with ", local.solver[i], ". Trying next algorithm \n", sep = "")
    } else {
      conv <- TRUE
      cat(" Solution converged with", local.solver[i], "\n")
      ss1 <- nlopt.ss$par
    }
    i <- i+1
  }

  if(nlopt.ss$convergence < 0 | all(nlopt.ss$par == ss0) | any(nlopt.ss$par <= 0)) {
    message("Consider relaxing one of the fixed sample size in the form c(x-.444, x+.444)")
    stop("Solution is not feasible. Change default settings", call.=FALSE)
  }

  col.names <- c(c(ifelse(!is.null(n1), ifelse(length(n1) == 2, "<n1<", "[n1]"),"n1"),
                   ifelse(!is.null(n2), ifelse(length(n2) == 2, "<n2<", "[n2]"),"n2"),
                   ifelse(!is.null(n3), ifelse(length(n3) == 2, "<n3<", "[n3]"),"n3"),
                   ifelse(!is.null(n4), ifelse(length(n4) == 2, "<n4<", "[n4]"),"n4"))[1:nlevels],
                 ifelse(!is.null(p), ifelse(length(p) == 2, "<p<", "[p]"),"p"),
                 ifelse(constrain == "cost","[cost]","cost"),
                 ifelse(constrain == "es","[es]","es"),
                 paste0(100 * round((1 - alpha), 2), "%lcl"),
                 paste0(100 * round((1 - alpha), 2), "%ucl"),
                 ifelse(constrain == "power","[power]","power"))

  #ss1 <- nlopt.ss$par
  # exact solution
  est.mlu.power <- .mlu.pwr(ss1)
  exact.cosa <- cbind(t(ss1), .min.cost(ss1), est.mlu.power[1], est.mlu.power[2], est.mlu.power[3], est.mlu.power[4])
  colnames(exact.cosa) <-  col.names

  # top-down rounded solution
  ss1[nlevels] <- round(ss1[nlevels])
  if(nlevels==2){
    ss1[1] <- round(prod(ss1[1:2]))/ss1[2]
  }else if(nlevels==3){
    ss1[2] <- round(prod(ss1[2:3]))/ss1[3]
    ss1[1] <- round(prod(ss1[1:3]))/prod(ss1[2:3])
  }else if(nlevels==4){
    ss1[3] <- round(prod(ss1[3:4]))/ss1[4]
    ss1[2] <- round(prod(ss1[2:4]))/prod(ss1[3:4])
    ss1[1] <- round(prod(ss1[1:4]))/prod(ss1[2:4])
  }
  if(nlevels==rlevel){
    ss1[nlevels+1] <- round(ss1[nlevels+1]*ss1[rlevel]) / ss1[rlevel]
  }else{
    ss1[nlevels+1] <- round(ss1[nlevels+1]*round(prod(ss1[rlevel:nlevels]))) / round(prod(ss1[rlevel:nlevels]))
  }

  est.mlu.power <- .mlu.pwr(ss1)
  round.cosa <- cbind(t(ss1), .min.cost(ss1), est.mlu.power[1], est.mlu.power[2], est.mlu.power[3], est.mlu.power[4])
  colnames(round.cosa) <-  col.names

  cat("\n Rounded solution: \n")
  cat(" ----------------------------------------------- \n")
  print(as.data.frame(round(round.cosa,3)), row.names=FALSE)
  cat(" ----------------------------------------------- \n")
  cat(" Per unit marginal costs: \n")
  cat(" Level 1 treatment:", cn1[1] ,
      "\n Level 1 control:", cn1[2], "\n")
  if(nlevels >= 2 & rlevel >= 2){
    cat(" Level 2 treatment:", cn2[1] ,
        "\n Level 2 control:", cn2[2], "\n")
  }else if(nlevels >= 2 & rlevel < 2){
    cat(" Level 2:", cn2[1], "\n")
  }
  if(nlevels >= 3 & rlevel >= 3){
    cat(" Level 3 treatment:", cn3[1] ,
        "\n Level 3 control:", cn3[2], "\n")
  }else if(nlevels >= 3 & rlevel < 3){
    cat(" Level 3:", cn3[1], "\n")
  }
  if(nlevels >= 4 & rlevel >= 4){
    cat(" Level 4 treatment:", cn4[1] ,
        "\n Level 4 control:", cn4[2], "\n")
  }else if(nlevels >= 4 & rlevel < 4){
    cat(" Level 4:", cn4[1], "\n")
  }
  cat(" ----------------------------------------------- \n")
  cat(" MDES = ", round(est.mlu.power[1], 3), " (with power = ", round(power,3)*100,
      ") \n power = ", round(est.mlu.power[4], 3), " (for ES = ", round(es,3), ") \n", sep="")
  cat(" ----------------------------------------------- \n")
  cat(" []: point constrained (fixed) \n <<: bound constrained \n")

  # output
  cosa.out <- list(parms = get("parms", parent.frame()),
                   exact.cosa = exact.cosa,
                   round.cosa = round.cosa
  )

  return(invisible(cosa.out))
}

# summarize mdes output
.summarize.mdes <- function(power, alpha, sse, df, two.tailed, mdes) {
  cat("\n Minimum detectable effect size: \n")
  cat(" --------------------------------------- \n")
  cat(" MDES is ", round(mdes[1], 3), " ", 100 * (1 - round(alpha, 2)),
      "% CI [", round(mdes[2], 3), ",", round(mdes[3], 3), "] with ",
      round(power,3)*100, "% power \n", sep="")
  cat(" --------------------------------------- \n")
  cat(" Degrees of freedom:", df, "\n Standardized standard error:", round(sse, 3),
      "\n Type I error rate:", alpha, "\n Type II error rate:", round(1 - power, 3),
      "\n Two-tailed test:", two.tailed)
}

# summarize power output
.summarize.power <- function(es, alpha, sse, df, two.tailed, power) {
  mlu <- .mdes(power, alpha, sse, df, two.tailed)
  cat("\n Statistical power: \n")
  cat(" --------------------------------------- \n")
  cat(" ", round(power, 3) * 100, "% power to detect an ES of ",
      round(mlu[1], 3), " ",  100 * (1 - round(alpha, 2)),
      "% CI [", round(mlu[2], 3), ",", round(mlu[3],3),"] \n", sep="")
  cat(" --------------------------------------- \n")
  cat(" Degrees of freedom:", df, "\n Standardized standard error:", round(sse, 3),
      "\n Type I error rate:", alpha, "\n Type II error rate:", round(1 - power, 3),
      "\n Two-tailed test:", two.tailed)
}
