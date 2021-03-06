# function to vectorize COSA solutions
vectorize.cosa <- function(design, args.grid, args.names = NULL, ordered = TRUE, ncase = 10L) {

  args.grid <- as.matrix(args.grid)
  n.cases <- nrow(args.grid)
  n.args <- ncol(args.grid)
  n.cosa <- length(design$cosa)

  if(is.null(args.names)) {
    if(isTRUE("Var" %in% substr(colnames(args.grid), start = 1, stop = 3))) {
      stop("Invalid argument names", call. = FALSE)
    } else {
      args.names <- colnames(args.grid)
    }
  }

  if(inherits(design, "cosa")){

    fun <- paste(class(design), collapse = ".")
    parms <- design$parms

    out <- matrix(ncol = n.cosa + n.args, nrow = n.cases)

    if(!isTRUE(ncol(args.grid) == n.args)){
      stop("'args.names' length is not consistent with number of columns in 'args.grid'", call. = FALSE)
    }

    if(isTRUE("cn" %in% substr(colnames(args.grid), start = 1, stop = 2))) {
      stop("Marginal costs cannot be vectorized", call. = FALSE)
    }

    if(isTRUE(parms[["round"]])){
      message("To reduce the vectorization time specify 'round = FALSE' in the main function")
    }

    if(n.cases > 50) {
      # to estimate lapsed time
      t1 <- Sys.time()
      capture.output(do.call(fun, parms))
      lapsed <- Sys.time() - t1
      cat("Estimated time for", n.cases, "runs:",
          round((n.cases * lapsed) / 60, 2) , "minutes \n")
    }

    if(all(!is.na(match(args.names, names(parms))))){

      for(i in 1:n.cases){
        parms[args.names] <- out[i, 1:n.args] <- args.grid[i,]
        capture.output({
          temp.out <- try(do.call(fun, parms)$cosa)
          if(inherits(temp.out, "try-error")) {
            out[i, (n.args + 1):(n.cosa + n.args)] <- NaN
          } else {
            out[i, (n.args + 1):(n.cosa + n.args)] <- temp.out
          }
        })
      }

      out <- cbind(1:n.cases, out)
      colnames(out) <- c("case", args.names, colnames(design$cosa))

      if(isTRUE(ordered)){
        out <- switch(parms$constrain,
                       "cost" = {out[order(out[, ncol(out)], decreasing = FALSE),]},
                       "power" = {out[order(out[, ncol(out) - 4], decreasing = TRUE),]},
                       "es" = {out[order(out[,  ncol(out) - 4], decreasing = TRUE),]}
                      )
        switch(parms$constrain,
               "cost" = cat("Top", ncase ,"solutions that have worst power rates \n"),
               "power" = cat("Top", ncase ,"solutions that have highest total cost \n"),
               "es" = cat("Top", ncase ,"solutions that have highest total cost \n"))
      }

      if(isFALSE(ordered)){
        ncase <- nrow(out)
        message("'ncase' argument is ignored because 'ordered = FALSE'")
        }

      if(isTRUE(ncase > nrow(out))){
        message("Number of cases to be printed exceeds arguments' grid length")
        ncase <- nrow(out)
      }

      return(round(out, 3)[1:ncase,])

    } else {

      stop("Invalid parameter(s)", call. = FALSE)

    }

  } else {

    stop("Only objects from COSA functions can be vectorized", call. = FALSE)

  }

}
