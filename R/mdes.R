mdes <- function(design = "ird1r1", power = .80, alpha = .05, two.tailed = TRUE,
                 rhots = NULL, k = c(-6, 6), dists = "normal",
                 rho = NULL, omega = NULL, r2t = NULL,
                 r2 = 0, n = 250, g = 0, p = .50) {

  message("\nExperimental function \n")
  if(is.null(design) || !is.character(design)) {
    stop("Incorrect value for argument 'design'", call. = FALSE)
  }
  block <- substr(design, nchar(design) - 1, nchar(design) - 1)
  rlevel <- as.numeric(substr(design, nchar(design), nchar(design)))
  nlevels <- as.numeric(substr(design, nchar(design) - 2, nchar(design) - 2))

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
  if (is.list(n)) {
    n <- unlist(n)
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

  cat("Function: mdes.", design, "() \n", sep = "")
  mdes.out <- switch(tolower(design),
         "ird1r1" = {
           if (n == 250) {
             cat(" Default sample size is 250 \n")
           }
           mdes.ird1r1(power = power, alpha = alpha, two.tailed = two.tailed,
                       rhots = rhots, k1 = k[1], k2 = k[2], dists = dists,
                       r21 = r2, g1 = g, p = p, n1 = n)
         },
         "bird2f1" = {
           mdes.bird2f1(power = power, alpha = alpha, two.tailed = two.tailed,
                       rhots = rhots, k1 = k[1], k2 = k[2], dists = dists,
                       r21 = r2, g1 = g, p = p, n1 = n[1], n2 = n[2])
         },
         "bird2r1" = {
           mdes.bird2r1(power = power, alpha = alpha, two.tailed = two.tailed,
                        rhots = rhots, k1 = k[1], k2 = k[2], dists = dists,
                        rho2 = rho, omega2 = omega, r21 = r2, r2t2 = r2t,
                        g2 = g, p = p, n1 = n[1], n2 = n[2])
         },
         "bird3r1" = {
           mdes.bird3r1(power = power, alpha = alpha, two.tailed = two.tailed,
                        rhots = rhots, k1 = k[1], k2 = k[2], dists = dists,
                        rho2 = rho[1], rho3 = rho[2], omega2 = omega[1], omega3 = omega[2],
                        r21 = r2, r2t2 = r2t[1], r2t3 = r2t[2], g3 = g,
                        p = p, n1 = n[1], n2 = n[2], n3 = n[3])
         },
         "bird4r1" = {
           mdes.bird4r1(power = power, alpha = alpha, two.tailed = two.tailed,
                        rhots = rhots, k1 = k[1], k2 = k[2], dists = dists,
                        rho2 = rho[1], rho3 = rho[2], rho4 = rho[3],
                        omega2 = omega[1], omega3 = omega[2], omega4 = omega[3],
                        r21 = r2, r2t2 = r2t[1], r2t3 = r2t[2], r2t4 = r2t[3],
                        g4 = g, p = p, n1 = n[1], n2 = n[2], n3 = n[3], n4 = n[4])
         },
         "crd2r2" = {
           mdes.crd2r2(power = power, alpha = alpha, two.tailed = two.tailed,
                        rhots = rhots, k1 = k[1], k2 = k[2], dists = dists,
                        rho2 = rho, r21 = r2[1], r22 = r2[2], g2 = g,
                        p = p, n1 = n[1], n2 = n[2])
         },
         "bcrd3f2" = {
           mdes.bcrd3f2(power = power, alpha = alpha, two.tailed = two.tailed,
                       rhots = rhots, k1 = k[1], k2 = k[2], dists = dists,
                       rho2 = rho, r21 = r2[1], r22 = r2[2], g2 = g,
                       p = p, n1 = n[1], n2 = n[2], n3 = n[3])
         },
         "bcrd3r2" = {
           mdes.bcrd3r2(power = power, alpha = alpha, two.tailed = two.tailed,
                        rhots = rhots, k1 = k[1], k2 = k[2], dists = dists,
                        rho2 = rho[1], rho3 = rho[2], omega3 = omega,
                        r21 = r2[1], r22 = r2[2], r2t3 = r2t,
                        g3 = g, p = p, n1 = n[1], n2 = n[2], n3 = n[3])
         },
         "bcrd4r2" = {
           mdes.bcrd4r2(power = power, alpha = alpha, two.tailed = two.tailed,
                        rhots = rhots, k1 = k[1], k2 = k[2], dists = dists,
                        rho2 = rho[1], rho3 = rho[2], rho4 = rho[3],
                        omega3 = omega[1], omega4 = omega[2],
                        r21 = r2[1], r22 = r2[2], r2t3 = r2t[1], r2t4 = r2t[2],
                        g4 = g, p = p, n1 = n[1], n2 = n[2], n3 = n[3], n4 = n[4])
         },
         "crd3r3" = {
           mdes.crd3r3(power = power, alpha = alpha, two.tailed = two.tailed,
                       rhots = rhots, k1 = k[1], k2 = k[2], dists = dists,
                       rho2 = rho[1], rho3 = rho[2], r21 = r2[1], r22 = r2[2], r23 = r2[3],
                       g3 = g, p = p, n1 = n[1], n2 = n[2], n3 = n[3])
         },
         "bcrd4f3" = {
           mdes.bcrd4f3(power = power, alpha = alpha, two.tailed = two.tailed,
                       rhots = rhots, k1 = k[1], k2 = k[2], dists = dists,
                       rho2 = rho[1], rho3 = rho[2], r21 = r2[1], r22 = r2[2], r23 = r2[3],
                       g3 = g, p = p, n1 = n[1], n2 = n[2], n3 = n[3], n4 = n[4])
         },
         "bcrd4r3" = {
           mdes.bcrd4r3(power = power, alpha = alpha, two.tailed = two.tailed,
                        rhots = rhots, k1 = k[1], k2 = k[2], dists = dists,
                        rho2 = rho[1], rho3 = rho[2], rho4 = rho[3], omega4 = omega,
                        r21 = r2[1], r22 = r2[2], r23 = r2[3], r2t4 = r2t,
                        g4 = g, p = p, n1 = n[1], n2 = n[2], n3 = n[3], n4 = n[4])
         },
         "crd4r4" = {
           mdes.crd4r4(power = power, alpha = alpha, two.tailed = two.tailed,
                       rhots = rhots, k1 = k[1], k2 = k[2], dists = dists,
                       rho2 = rho[1], rho3 = rho[2], rho4 = rho[3],
                       r21 = r2[1], r22 = r2[2], r23 = r2[3], r24 = r2[4],
                       g4 = g, p = p, n1 = n[1], n2 = n[2], n3 = n[3], n4 = n[4])
         },
         stop("Incorrect design. Use one of the: 'ird1r1', 'bird2f1', 'bird2r1', 'bird3r1', 'bird4r1',
              'crd2r2', 'bcrd3f2', 'bcrd3r2', 'bcrd4r2', 'crd3r3', 'bcrd4f3', 'bcrd4r3', 'crd4r4'", call. = FALSE)
         )

  return(invisible(mdes.out))
}



# examples

# mdes.crd4r4(rho4 = .05, rho3 = .05, rho2 = .10, n1 = 10, n2 = 2, n3 = 3, n4 = 20)
# mdes(design = "crd4r4", rho = c(.10, .05, .05), n = c(10, 2, 3, 20), r2 = c(0, 0, 0, 0))

# mdes.crd3r3(rho3 = .06, rho2 = .17, n1 = 15, n2 = 3, n3 = 60)
# mdes(design = "crd3r3", rho = c(.17, .06), n = c(15, 3, 60), r2 = c(0, 0, 0))

# mdes.crd2r2(rho2=.20, n1 = 4, n2 = 20)
# mdes(design = "crd2r2", rho = .20, n = c(4, 20), r2 = c(0, 0))

# mdes.ird1r1(n = 400)
# mdes(design = "ird1r1", n = 400, r2 = 0)
