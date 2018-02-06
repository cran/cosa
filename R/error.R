# error handler
.error.handler <- function(x){

  # exclude default NULL arguments
  names.x <- names(x)
  idx.notnull <- match(names(lapply(x, is.null)[!lapply(x, is.null) == TRUE]), names.x)
  parms.notnull <- x[idx.notnull]

  # redefine the check list
  names.x <- names(parms.notnull)
  x <-  parms.notnull

  # validity check for sample sizes
  idx.N <- match(c("n1","n2","n3","n4"),  names.x)
  if(any(unlist(x[idx.N[!is.na(idx.N) &
                        lapply(x[idx.N], is.null) == FALSE]]) <= 0)){
    stop("Negative sample size", call. = FALSE)
  }

  # validity check for number of covariates
  idx.g <- match(c("g", "g1", "g2", "g3", "g4"),  names.x)
  if(any(x[idx.g[!is.na(idx.g)]] < 0)){
    stop("Negative number of covariates", call. = FALSE)
  }

  # validity check for variance parameters, proportions, and probabilities
  idx.var <- match(c("r21","r22","r23", "r24", "r2t2", "r2t3", "r2t4",
                     "rho2", "rho3", "rho4", "omega2", "omega3", "omega4",
                     "alpha", "power"),  names.x)
  if(any(x[idx.var[!is.na(idx.var)]] < 0) | any(x[idx.var[!is.na(idx.var)]] > 1)){
    stop("Argument out of range when it should be between 0 and 1", call. = FALSE)
  }
  if(!is.null(x$p) & any(x$p < .01 | x$p > .99 )){
    stop("p should be between .01 and .99", call. = FALSE)
  }
  idx.r2 <- match(c("r2", "r21","r22","r23", "r24", "r2t2", "r2t3", "r2t4"),  names.x)
  if(any(x[idx.r2[!is.na(idx.r2)]] > 0) & any(x[idx.g[!is.na(idx.g)]] == 0)){
    err.r2 <- names.x[idx.r2[!is.na(idx.r2)]][x[idx.r2[!is.na(idx.r2)]] > 0]
    err.g <- names.x[idx.g[!is.na(idx.g)]][x[idx.g[!is.na(idx.g)]] == 0]
    if(any(substr(err.r2, 4, 4) == substr(err.g, 2, 2))){
      warning("R-squared value for a level may not be greater than zero
              when the number of covariates at that level [other than blocking variables] is zero", call. = FALSE)
    }
  }else if(any(x[idx.r2[!is.na(idx.r2)]] == 0) & any(x[idx.g[!is.na(idx.g)]] > 0 )){
    err.r2 <- names.x[idx.r2[!is.na(idx.r2)]][x[idx.r2[!is.na(idx.r2)]] == 0]
    err.g <- names.x[idx.g[!is.na(idx.g)]][x[idx.g[!is.na(idx.g)]] > 0]
    if(any(substr(err.r2, 4, 4) == substr(err.g, 2, 2))){
      warning("R-squared value for a level may not be zero
              when the number of covariates at that level [other than blocking variables] is greater than zero ", call. = FALSE)
    }
  }

  # warn for negative treatment effect
  if(any(x$es < 0)){
    warning("Negative ES", call. = FALSE)
  }

  # validty check for two-tailed test
  if(x$two.tailed != TRUE & x$two.tailed != FALSE){
    stop("Non-logical value for two.tailed", call. = FALSE)
  }

}#.error.handler()


