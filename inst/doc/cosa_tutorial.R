## ----setup, include=FALSE, warning=FALSE---------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(cosa)

## ----eval = FALSE--------------------------------------------------------
#  install.packages("cosa")
#  library(cosa)

## ---- message=FALSE, fig.width=7, fig.height=5, warning=FALSE------------
# cost constrained - optimize p and n2
# CRT('order = 0' or 'rhots = 0')
crt <- cosa.crd2(order = 0,
                 constrain = "cost", cost = 12500, 
                 cn1 = c(5, 2), cn2 = c(50, 20),
                 es = .20,  power = .80, rho2 = .20,
                 g2 = 1, r21 = .20, r22 = .30, 
                 p = NULL, n1 = 20,  n2 = NULL)
# comparisons to CRDs
# CRD w/ linear score variable interacting with the treatment
crd1 <- cosa.crd2(order = 1, interaction = TRUE,
                  constrain = "cost", cost = 12500,
                  cn1 = c(5, 2), cn2 = c(50, 20),
                  es = .20,  power = .80, rho2 = .20,
                  g2 = 1, r21 = .20, r22 = .30,
                  p = .386, n1 = 24,  n2 = NULL)
# CRD w/ quadratic score variable interacting with the treatment
crd2 <- cosa.crd2(order = 2, interaction = TRUE,
                  constrain = "cost", cost = 12500,
                  cn1 = c(5, 2), cn2 = c(50, 20),
                  es = .20,  power = .80, rho2 = .20,
                  g2 = 1, r21 = .20, r22 = .30,
                  p = .386, n1 = 24,  n2 = NULL)

# example plots
par(mfrow = c(2, 3), mai = c(.6, .6, .6, .2))
# compare minimum detectable effect size and 95% CI
plot(crt, ypar = "mdes", xpar = "n2",
     ylim = c(.10, .90), xlim = c(10, 800), 
     ylab = "MDES (with Power = .80)", xlab = "Number of Clusters",
     main = expression(CRT), locate = TRUE)
plot(crd1, ypar = "mdes", xpar = "n2",
     ylim = c(.10, .90), xlim = c(10, 800),
     ylab = "MDES (with Power = .80)", xlab = "Number of Clusters",
     main = expression(CRD (S + TS)), locate = TRUE)
plot(crd2, ypar = "mdes", xpar = "n2",
     ylim = c(.10, .90), xlim = c(10, 800),
     ylab = "MDES (with Power = .80)", xlab = "Number of Clusters",
     main = expression(CRD (S + S^2 + TS + TS^2)), locate = TRUE)
# compare statistical power 
plot(crt, ypar = "power", xpar = "n2",
     ylim = c(.10, .85), xlim = c(10, 800), 
     ylab = "Power (for ES = .20)", xlab = "Number of Clusters",
     main = expression(CRT), locate = TRUE)
plot(crd1, ypar = "power", xpar = "n2",
     ylim = c(.10, .85), xlim = c(10, 800), 
     ylab = "Power (for ES = .20)", xlab = "Number of Clusters",
     main = expression(CRD(S + TS)), locate = TRUE)
plot(crd2, ypar = "power", xpar = "n2",
     ylim = c(.10, .85), xlim = c(10, 800), 
     ylab = "Power (for ES = .20)", xlab = "Number of Clusters",
     main = expression(CRD (S + S^2 + TS + TS^2)), locate = TRUE)


## ---- message=FALSE, warning=FALSE---------------------------------------
# cost constrained - optimize p and n2
# CRT('order = 0' or 'rhots = 0')
cosa.crd2(order = 0,
          constrain = "es", es = .20,  
          cn1 = c(5, 2), cn2 = c(50, 20),
          power = .80, rho2 = .20,
          g2 = 1, r21 = .20, r22 = .30, 
          p = NULL, n1 = 20,  n2 = NULL)
# comparisons to CRDs
# CRD w/ linear score variable interacting with the treatment
cosa.crd2(order = 1, interaction = TRUE,
          constrain = "es", es = .20,  
          cn1 = c(5, 2), cn2 = c(50, 20),
          power = .80, rho2 = .20,
          g2 = 1, r21 = .20, r22 = .30,
          p = .389, n1 = 24,  n2 = NULL)
# CRD w/ quadratic score variable interacting with the treatment
cosa.crd2(order = 2, interaction = TRUE,
          constrain = "es", es = .20,  
          cn1 = c(5, 2), cn2 = c(50, 20),
          power = .80, rho2 = .20,
          g2 = 1, r21 = .20, r22 = .30,
          p = .389, n1 = 24,  n2 = NULL)


## ---- message=FALSE, warning=FALSE---------------------------------------
# cost constrained - optimize p and n2
# CRT('order = 0' or 'rhots = 0')
cosa.crd2(order = 0,
          constrain = "power", power = .80, 
          cn1 = c(5, 2), cn2 = c(50, 20),
          es = .20, rho2 = .20,
          g2 = 1, r21 = .20, r22 = .30, 
          p = NULL, n1 = 20,  n2 = NULL)
# comparisons to CRDs
# CRD w/ linear score variable interacting with the treatment
cosa.crd2(order = 1, interaction = TRUE,
          constrain = "power", power = .80, 
          cn1 = c(5, 2), cn2 = c(50, 20),
          es = .20, rho2 = .20,
          g2 = 1, r21 = .20, r22 = .30,
          p = .389, n1 = 24,  n2 = NULL)
# CRD w/ quadratic score variable interacting with the treatment
cosa.crd2(order = 2, interaction = TRUE,
          constrain = "power", power = .80,  
          cn1 = c(5, 2), cn2 = c(50, 20),
          es = .20, rho2 = .20,
          g2 = 1, r21 = .20, r22 = .30,
          p = .389, n1 = 24,  n2 = NULL)


