cosa <- function(design = "ird1r1",
                 cn = as.list(rep(0, length(n))), cost = NULL,
                 n = list(NULL), p = NULL, n0 = rep(5 + g, length(n)), p0 = .499,
                 constrain = "power", local.solver = c("LBFGS", "SLSQP", "MMA", "COBYLA"),
                 power = .80, es = .25, alpha = .05, two.tailed = TRUE,
                 rhots = NULL, k = c(-6, 6), dists = "normal",
                 rho = NULL, omega = NULL, r2t = NULL,
                 r2 = 0,  g = 0) {

  message("\nExperimental function \n")
  if(is.null(design) || !is.character(design)) {
    stop("Incorrect value for argument 'design'", call. = FALSE)
  }
  block <- substr(design, nchar(design) - 1, nchar(design) - 1)
  rlevel <- as.numeric(substr(design, nchar(design), nchar(design)))
  nlevels <- as.numeric(substr(design, nchar(design) - 2, nchar(design) - 2))

  if (!is.list(cn)) {
    stop("'cn' must be a list", call. = FALSE)
  }
  if (!is.list(n)) {
    stop("'n' must be a list", call. = FALSE)
  }
  if (is.list(local.solver)) {
    local.solver <- unlist(local.solver)
  }
  if (is.list(n0)) {
    n0 <- unlist(n0)
  }
  if (is.list(rho)) {
    rho <- unlist(rho)
  }
  if (is.list(omega)) {
    omega <- unlist(omega)
  }
  if (is.list(r2)) {
    r2 <- unlist(r2)
  }
  if (is.list(r2t)) {
    r2t <- unlist(r2t)
  }
  if (is.list(k)) {
    k <- unlist(k)
  }

  if (any(c(length(omega), length(r2t)) != nlevels - rlevel)) {
    if(block == "r") {
      stop("Length of 'omega' or 'r2t' is not consistent with 'design'", call. = FALSE)
    }
  }
  if (length(r2) != rlevel) {
    stop("Length of 'r2' is not consistent with 'design'", call. = FALSE)
  }
  if (length(n) != nlevels) {
    stop("Length of 'n' is not consistent with 'design'", call. = FALSE)
  }
  if (length(rho) != nlevels - 1) {
    if(block == "r") {
      stop("Length of 'rho' is not consistent with 'design'", call. = FALSE)
    }
  }

  cat("Function: cosa.", design, "() \n", sep = "")
  cosa.out <- switch(tolower(design),
         "ird1r1" = {
           cosa.ird1r1(n0 = n0, p0 = p0,
                       cn1 = cn[[1]], cost = cost,
                       constrain = constrain, local.solver = local.solver,
                       es = es, alpha = alpha, two.tailed = two.tailed,
                       rhots = rhots, k1 = k[1], k2 = k[2], dists = dists,
                       r21 = r2, g1 = g, p = p, n1 = n[[1]])
         },
         "bird2f1" = {
           cosa.bird2f1(n0 = n0, p0 = p0,
                       cn1 = cn[[1]], cn2 = cn[[2]], cost = cost,
                       constrain = constrain, local.solver = local.solver,
                       es = es, alpha = alpha, two.tailed = two.tailed,
                       rhots = rhots, k1 = k[1], k2 = k[2], dists = dists,
                       r21 = r2, g1 = g, p = p, n1 = n[[1]], n2 = n[[2]])
         },
         "bird2r1" = {
           cosa.bird2r1(n0 = n0, p0 = p0,
                        cn1 = cn[[1]], cn2 = cn[[2]], cost = cost,
                        constrain = constrain, local.solver = local.solver,
                        es = es, alpha = alpha, two.tailed = two.tailed,
                        rhots = rhots, k1 = k[1], k2 = k[2], dists = dists,
                        rho2 = rho, omega2 = omega, r21 = r2, r2t2 = r2t,
                        g2 = g, p = p, n1 = n[[1]], n2 = n[[2]])
         },
         "bird3r1" = {
           cosa.bird3r1(n0 = n0, p0 = p0,
                        cn1 = cn[[1]], cn2 = cn[[2]], cn3 = cn[[3]], cost = cost,
                        constrain = constrain, local.solver = local.solver,
                        es = es, alpha = alpha, two.tailed = two.tailed,
                        rhots = rhots, k1 = k[1], k2 = k[2], dists = dists,
                        rho2 = rho[1], rho3 = rho[2], omega2 = omega[1], omega3 = omega[2],
                        r21 = r2, r2t2 = r2t[1], r2t3 = r2t[2], g3 = g,
                        p = p, n1 = n[[1]], n2 = n[[2]], n3 = n[[3]])
         },
         "bird4r1" = {
           cosa.bird4r1(n0 = n0, p0 = p0,
                        cn1 = cn[[1]], cn2 = cn[[2]], cn3 = cn[[3]], cn4 = cn[[4]], cost = cost,
                        constrain = constrain, local.solver = local.solver,
                        es = es, alpha = alpha, two.tailed = two.tailed,
                        rhots = rhots, k1 = k[1], k2 = k[2], dists = dists,
                        rho2 = rho[1], rho3 = rho[2], rho4 = rho[3],
                        omega2 = omega[1], omega3 = omega[2], omega4 = omega[3],
                        r21 = r2, r2t2 = r2t[1], r2t3 = r2t[2], r2t4 = r2t[3],
                        g4 = g, p = p, n1 = n[[1]], n2 = n[[2]], n3 = n[[3]], n4 = n[[4]])
         },
         "crd2r2" = {
           cosa.crd2r2(n0 = n0, p0 = p0,
                       cn1 = cn[[1]], cn2 = cn[[2]], cost = cost,
                       constrain = constrain, local.solver = local.solver,
                       es = es, alpha = alpha, two.tailed = two.tailed,
                       rhots = rhots, k1 = k[1], k2 = k[2], dists = dists,
                       rho2 = rho, r21 = r2[1], r22 = r2[2], g2 = g,
                       p = p, n1 = n[[1]], n2 = n[[2]])
         },
         "bcrd3f2" = {
           cosa.bcrd3f2(n0 = n0, p0 = p0,
                       cn1 = cn[[1]], cn2 = cn[[2]], cn3 = cn[[3]], cost = cost,
                       constrain = constrain, local.solver = local.solver,
                       es = es, alpha = alpha, two.tailed = two.tailed,
                       rhots = rhots, k1 = k[1], k2 = k[2], dists = dists,
                       rho2 = rho, r21 = r2[1], r22 = r2[2], g2 = g,
                       p = p, n1 = n[[1]], n2 = n[[2]], n3 = n[[3]])
         },
         "bcrd3r2" = {
           cosa.bcrd3r2(n0 = n0, p0 = p0,
                        cn1 = cn[[1]], cn2 = cn[[2]], cn3 = cn[[3]], cost = cost,
                        constrain = constrain, local.solver = local.solver,
                        es = es, alpha = alpha, two.tailed = two.tailed,
                        rhots = rhots, k1 = k[1], k2 = k[2], dists = dists,
                        rho2 = rho[1], rho3 = rho[2], omega3 = omega,
                        r21 = r2[1], r22 = r2[2], r2t3 = r2t,
                        g3 = g, p = p, n1 = n[[1]], n2 = n[[2]], n3 = n[[3]])
         },
         "bcrd4r2" = {
           cosa.bcrd4r2(n0 = n0, p0 = p0,
                        cn1 = cn[[1]], cn2 = cn[[2]], cn3 = cn[[3]], cn4 = cn[[4]], cost = cost,
                        constrain = constrain, local.solver = local.solver,
                        es = es, alpha = alpha, two.tailed = two.tailed,
                        rhots = rhots, k1 = k[1], k2 = k[2], dists = dists,
                        rho2 = rho[1], rho3 = rho[2], rho4 = rho[3],
                        omega3 = omega[1], omega4 = omega[2],
                        r21 = r2[1], r22 = r2[2], r2t3 = r2t[1], r2t4 = r2t[2],
                        g4 = g, p = p, n1 = n[[1]], n2 = n[[2]], n3 = n[[3]], n4 = n[[4]])
         },
         "crd3r3" = {
           cosa.crd3r3(n0 = n0, p0 = p0,
                       cn1 = cn[[1]], cn2 = cn[[2]], cn3 = cn[[3]], cost = cost,
                       constrain = constrain, local.solver = local.solver,
                       es = es, alpha = alpha, two.tailed = two.tailed,
                       rhots = rhots, k1 = k[1], k2 = k[2], dists = dists,
                       rho2 = rho[1], rho3 = rho[2], r21 = r2[1], r22 = r2[2], r23 = r2[3],
                       g3 = g, p = p, n1 = n[[1]], n2 = n[[2]], n3 = n[[3]])
         },
         "bcrd4f3" = {
           cosa.bcrd4f3(n0 = n0, p0 = p0,
                       cn1 = cn[[1]], cn2 = cn[[2]], cn3 = cn[[3]], cn4 = cn[[4]], cost = cost,
                       constrain = constrain, local.solver = local.solver,
                       es = es, alpha = alpha, two.tailed = two.tailed,
                       rhots = rhots, k1 = k[1], k2 = k[2], dists = dists,
                       rho2 = rho[1], rho3 = rho[2], r21 = r2[1], r22 = r2[2], r23 = r2[3],
                       g3 = g, p = p, n1 = n[[1]], n2 = n[[2]], n3 = n[[3]], n4 = n[[4]])
         },
         "bcrd4r3" = {
           cosa.bcrd4r3(n0 = n0, p0 = p0,
                        cn1 = cn[[1]], cn2 = cn[[2]], cn3 = cn[[3]], cn4 = cn[[4]], cost = cost,
                        constrain = constrain, local.solver = local.solver,
                        es = es, alpha = alpha, two.tailed = two.tailed,
                        rhots = rhots, k1 = k[1], k2 = k[2], dists = dists,
                        rho2 = rho[1], rho3 = rho[2], rho4 = rho[3], omega4 = omega,
                        r21 = r2[1], r22 = r2[2], r23 = r2[3], r2t4 = r2t,
                        g4 = g, p = p, n1 = n[[1]], n2 = n[[2]], n3 = n[[3]], n4 = n[[4]])
         },
         "crd4r4" = {
           cosa.crd4r4(n0 = n0, p0 = p0,
                       cn1 = cn[[1]], cn2 = cn[[2]], cn3 = cn[[3]], cn4 = cn[[4]], cost = cost,
                       constrain = constrain, local.solver = local.solver,
                       es = es, alpha = alpha, two.tailed = two.tailed,
                       rhots = rhots, k1 = k[1], k2 = k[2], dists = dists,
                       rho2 = rho[1], rho3 = rho[2], rho4 = rho[3],
                       r21 = r2[1], r22 = r2[2], r23 = r2[3], r24 = r2[4],
                       g4 = g, p = p, n1 = n[[1]], n2 = n[[2]], n3 = n[[3]], n4 = n[[4]])
         },
         stop("Incorrect design. Use one of the: 'ird1r1', 'bird2f1', 'bird2r1', 'bird3r1', 'bird4r1',
              'crd2r2', 'bcrd3f2', 'bcrd3r2', 'bcrd4r2', 'crd3r3', 'bcrd4f3', 'bcrd4r3', 'crd4r4'", call. = FALSE)
         )

  return(invisible(cosa.out))
}

# examples
# cosa(design = "ird1r1", n = list(NULL), r2 = 0, p = .50)

# cosa.crd4r4(rho4 = .05, rho3 = .05, rho2 = .10, n1 = 10, n2 = 2, p =.50)
# cosa(design = "crd4r4", rho = c(.10, .05, .05), n = list(10, 2, NULL, NULL), r2 = c(0, 0, 0, 0), p = .50)

# cosa.crd4r4(rhots = 0, rho4 = .05, rho3 = .05, rho2 = .10, n1 = 10, n2 = 2, n3 = 20)
# cosa(design = "crd4r4", rhots = 0,  rho = c(.10, .05, .05), n = list(10, 2, 20, NULL), r2 = c(0, 0, 0, 0))

# cosa.bcrd4r2(rhots = 0, rho4 = .05, rho3 = .05, rho2 = .10,
#             omega3 = .20, omega4 = .20, n1 = 10, n2 = 2, n3 = c(10, 12))
# cosa(design = "bcrd4r2", rhots = 0,
#     rho = c(.10, .05, .05),
#     omega = c(.20, .20),
#     r2t = c(0, 0),
#     n = list(10, 2, c(10, 12), NULL),
#     r2 = c(0, 0))

# cosa.crd3r3(rho3 = .06, rho2 = .17, n1 = 15, p = .50)
# cosa(design = "crd3r3", rho = c(.17, .06), n = list(15, NULL, NULL), r2 = c(0, 0, 0), p =.50)
# cosa(design = "crd3r3", rhots = 0,
#      cn =list(c(5, 3), c(10,5), c(50,20)), cost = 25000, constrain = "cost",
#      rho = c(.17, .06), n = list(15, NULL, NULL), r2 = c(0, 0, 0), p =.50)

# cosa.crd3r3(rhots = 0,  rho3 = .06, rho2 = .17, n1 = 10, n2 = 3)
# cosa(design = "crd3r3", rhots = 0, rho = c(.17, .06), n = list(10, 3, NULL), r2 = c(0, 0, 0))

# cosa.crd2r2(rhots = 0, rho2=.20, n1 = 20)
# cosa(design = "crd2r2", rhots = 0, rho = .20, n = list(20, NULL), r2 = c(0, 0))

