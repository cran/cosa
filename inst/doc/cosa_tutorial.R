## ----setup, include=FALSE, warning=FALSE---------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(cosa)

## ----eval = FALSE--------------------------------------------------------
#  install.packages("cosa")
#  library(cosa)

## ---- message=FALSE, fig.width=7, fig.height=5---------------------------
# cost constrained - optimize p and n2
# CRT('order = 0' or 'rhots = 0')
crt <- cosa.crd2r2(order = 0,
                   constrain = "cost", cost = 12500, 
                   cn1 = c(5, 2), cn2 = c(50, 20),
                   es = .20,  power = .80, rho2 = .20,
                   g2 = 5, r21 = .20, r22 = .30, 
                   p = NULL, n1 = 20,  n2 = NULL)
# comparisons to CRDs
# CRD w/ linear score variable
crd1 <- cosa.crd2r2(order = 1,
                    constrain = "cost", cost = 12500,
                    cn1 = c(5, 2), cn2 = c(50, 20),
                    es = .20,  power = .80, rho2 = .20,
                    g2 = 5, r21 = .20, r22 = .30,
                    p = .386, n1 = 24,  n2 = NULL)
# CRD w/ linear + quadratic score variable
crd2 <- cosa.crd2r2(order = 2,
                    constrain = "cost", cost = 12500,
                    cn1 = c(5, 2), cn2 = c(50, 20),
                    es = .20,  power = .80, rho2 = .20,
                    g2 = 5, r21 = .20, r22 = .30,
                    p = .386, n1 = 24,  n2 = NULL)

# example plots
par(mfrow = c(2, 3), mai = c(.6, .6, .6, .2))
# compare minimum detectable effect size and 95% CI
plot(crt, ypar = "mdes", xpar = "n2",
     ylim = c(.10, .50), xlim = c(10, 200), 
     ylab = "MDES (with Power = .80)", xlab = "Number of Clusters",
     main = expression(CRT), locate = TRUE)
plot(crd1, ypar = "mdes", xpar = "n2",
     ylim = c(.10, .50), xlim = c(10, 200),
     ylab = "MDES (with Power = .80)", xlab = "Number of Clusters",
     main = expression(CRD(S)), locate = TRUE)
plot(crd2, ypar = "mdes", xpar = "n2",
     ylim = c(.10, .50), xlim = c(10, 200),
     ylab = "MDES (with Power = .80)", xlab = "Number of Clusters",
     main = expression(CRD (S + S^2)), locate = TRUE)
# compare statistical power 
plot(crt, ypar = "power", xpar = "n2",
     ylim = c(.10, .85), xlim = c(10, 200), 
     ylab = "Power (for ES = .20)", xlab = "Number of Clusters",
     main = expression(CRT), locate = TRUE)
plot(crd1, ypar = "power", xpar = "n2",
     ylim = c(.10, .85), xlim = c(10, 200), 
     ylab = "Power (for ES = .20)", xlab = "Number of Clusters",
     main = expression(CRD(S)), locate = TRUE)
plot(crd2, ypar = "power", xpar = "n2",
     ylim = c(.10, .85), xlim = c(10, 200), 
     ylab = "Power (for ES = .20)", xlab = "Number of Clusters",
     main = expression(CRD (S + S^2)), locate = TRUE)


## ---- message=FALSE------------------------------------------------------
# cost constrained - optimize p and n2
# CRT('order = 0' or 'rhots = 0')
cosa.crd2r2(order = 0,
            constrain = "es", es = .20,  
            cn1 = c(5, 2), cn2 = c(50, 20),
            power = .80, rho2 = .20,
            g2 = 5, r21 = .20, r22 = .30, 
            p = NULL, n1 = 20,  n2 = NULL)
# comparisons to CRDs
# CRD w/ linear score variable
cosa.crd2r2(order = 1,
            constrain = "es", es = .20,  
            cn1 = c(5, 2), cn2 = c(50, 20),
            power = .80, rho2 = .20,
            g2 = 5, r21 = .20, r22 = .30,
            p = .386, n1 = 24,  n2 = NULL)
# CRD w/ linear + quadratic score variable
cosa.crd2r2(order = 2,
            constrain = "es", es = .20,  
            cn1 = c(5, 2), cn2 = c(50, 20),
            power = .80, rho2 = .20,
            g2 = 5, r21 = .20, r22 = .30,
            p = .386, n1 = 24,  n2 = NULL)


## ---- message=FALSE------------------------------------------------------
# cost constrained - optimize p and n2
# CRT('order = 0' or 'rhots = 0')
cosa.crd2r2(order = 0,
            constrain = "power", power = .80, 
            cn1 = c(5, 2), cn2 = c(50, 20),
            es = .20, rho2 = .20,
            g2 = 5, r21 = .20, r22 = .30, 
            p = NULL, n1 = 20,  n2 = NULL)
# comparisons to CRDs
# CRD w/ linear score variable
cosa.crd2r2(order = 1,
            constrain = "power", power = .80, 
            cn1 = c(5, 2), cn2 = c(50, 20),
            es = .20, rho2 = .20,
            g2 = 5, r21 = .20, r22 = .30,
            p = .386, n1 = 24,  n2 = NULL)
# CRD w/ linear + quadratic score variable
cosa.crd2r2(order = 2,
            constrain = "power", power = .80,  
            cn1 = c(5, 2), cn2 = c(50, 20),
            es = .20, rho2 = .20,
            g2 = 5, r21 = .20, r22 = .30,
            p = .386, n1 = 24,  n2 = NULL)


## ------------------------------------------------------------------------
# fix 'p = .50' to inpsect cost for the balanced case
cosa.crd2r2(order = 0,
            constrain = "power", power = .80, 
            cn1 = c(5, 2), cn2 = c(50, 20),
            es = .20, rho2 = .20,
            g2 = 5, r21 = .20, r22 = .30, 
            p = .50, n1 = 20,  n2 = NULL)

