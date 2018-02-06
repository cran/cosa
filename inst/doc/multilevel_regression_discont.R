## ----setup, include=FALSE, warning=FALSE---------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(cosa)

## ---- message=FALSE------------------------------------------------------
# power constrained - optimize n1 and n2
cosa.bird2r1(constrain = "power", power = .80,
             es = .20, rho2 = .20,  omega2 = .20,
             r21 = .20, r2t2 = .20, g2 = 1, p = .50,
             n1 = c(10, 40))

# cost constrained - optimize n1 and n2
cosa.bird2r1(constrain = "cost", cost = 2000,
             cn1 = c(5, 4), cn2 = 10,
             es = .20, rho2 = .20, omega2 = .20,
             r21 = .20, r2t2 = .20, g2 = 1, p = .50,
             n1 = c(10, 40))

## ---- message=FALSE------------------------------------------------------
# power constrained - optimize n1 and n3
cosa.bird3r1(constrain = "power", power = .80,
             es = .20, rho2 = .20, rho3 = .10,  omega2 = .20, omega3 = .10,
             r21 = .20, r2t2 = .20, r2t3 = .20, g3 = 1, p = .50,
             n1 = c(10, 40), n2 = 3)

# cost constrained - optimize n1 and n3
cosa.bird3r1(constrain = "cost", cost = 5000,
             cn1 = c(5, 4), cn2 = 10, cn3 = 30,
             es = .20, rho2 = .20, rho3 = .10, omega2 = .20, omega3 = .10,
             r21 = .20, r2t2 = .20, r2t3 = .20, g3 = 1, p = .50,
             n1 = c(10, 40), n2 = 3)

## ---- message=FALSE------------------------------------------------------
# cost constrained - optimize n1, and n3 
cosa.bird4r1(constrain = "cost", cost = 50000,
             cn1 = c(5, 4), cn2 = 10, cn3 = 30, cn4 = 50,
             es = .20, rho2 = .20, rho3 = .10, rho4 = .05,
             omega2 = .20, omega3 = .10, omega4 = .10,
             r21 = .20, r2t2 = .20, r2t3 = .20, r2t4 = .20,
             g4 = 1, p = .50,
             n1 = c(10, 40), n2 = 3)


## ---- message=FALSE------------------------------------------------------
# cost constrained - optimize n1 and n2
cosa.crd2r2(constrain = "cost", cost = 50000,
            cn1 = 5, cn2 = c(10, 8),
            es = .20, rho2 = .20, r21 = .20, r22 = .30,
            g2 = 1, p = .50)

## ---- message=FALSE------------------------------------------------------
# cost constrained
cosa.bcrd3r2(constrain = "cost", cost = 50000,
             cn1 = c(5, 4), cn2 = 10, cn3 = 30,
             es = .20, rho2 = .20, rho3 = .10, omega3 = .10,
             r21 = .20, r22 = .20, r2t3 = .20,
             g3 = 1, p = .50, n1 = c(10, 40), n2 = 3)

## ---- message=FALSE------------------------------------------------------
# cost constrained
cosa.bcrd4r2(constrain = "cost", cost = 80000,
             cn1 = c(5, 4), cn2 = 10, cn3 = 30, cn4 = 50,
             es = .20, rho2 = .20, rho3 = .10, rho4 = .05,
             omega3 = .10, omega4 = .20,
             r21 = .20, r22 = .20, r2t3 = .20, r2t4 = .20,
             g4 = 1, p = .50,
             n1 = 10, n2 = 3)

## ---- message=FALSE------------------------------------------------------
# cost constrained - optimize n1 and n3
cosa.crd3r3(constrain = "cost", cost = 100000,
            cn1 = 5, cn2 = c(10, 5), cn3 = c(30, 10),
            es = .20, rho2 = .20, rho3 = .10,
            r21 = .20, r22 = .30, r23 = .40,
            g3 = 1, p = .50,
            n1 = c(10, 40), n2 = 3)

## ---- message=FALSE------------------------------------------------------
# cost constrained - optimize n1, n3 and n4
cosa.bcrd4r3(constrain = "cost", cost = 100000,
             cn1 = c(5, 4), cn2 = 10, cn3 = 30, cn4 = 50,
             es = .20, rho2 = .20, rho3 = .10, rho4 = .05, omega4 = .20,
             r21 = .20, r22 = .20, r23 = .20, r2t4 = .20,
             g4 = 1, p = .50,
             n1 = 25, n2 = 3)

## ---- message=FALSE------------------------------------------------------
# cost constrained - optimize n1, n3 and n4
cosa.crd4r4(constrain = "cost", cost = 100000,
            cn1 = 5, cn2 = c(10, 5), cn3 = c(30, 10), cn4 = c(50,30),
            es = .20, rho2 = .20, rho3 = .10, rho4 = .05,
            r21 = .20, r22 = .30, r23 = .40, r24 = .50,
            g4 = 1, p = .50,
            n1 = 10, n2 = 3)

