## ----setup, include=FALSE, warning=FALSE---------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(cosa)

## ---- message=FALSE------------------------------------------------------
# cost constrained - optimize p and n2
# Specifying rhots as 0 allows COSA in CRT
crt <- cosa.crd2r2(rhots = 0,
                   constrain = "cost", cost = 15000,
                   cn1 = c(5, 2), cn2 = c(50, 20),
                   es = .20, rho2 = .20,
                   g2 = 5, r21 = .20, r22 = .30,
                   p = NULL, n1 = 24,  n2 = NULL)

## ---- message=FALSE------------------------------------------------------
# cost constrained - optimize n2
crd <- cosa.crd2r2(constrain = "cost", cost = 15000,
                   cn1 = c(5, 2), cn2 = c(50, 20),
                   es = .20, rho2 = .20,
                   g2 = 5, r21 = .20, r22 = .30,
                   p = .387, n1 = 24,  n2 = NULL)

## ---- message=FALSE, fig.width=7, fig.height=5---------------------------
par(mfrow = c(1, 2), mai = c(1, .9, .5, .2))

plot(crt, ypar = "mdes", xpar = "n2",
     ylim = c(0, 1), xlim = c(10, 200), 
     ylab = "MDES (Power = .80)", xlab = "Number of Clusters",
     main = "CRT", locate = TRUE)
plot(crd, ypar = "mdes", xpar = "n2",
     ylim = c(0, 1), xlim = c(10, 200),
     ylab = "MDES (Power = .80)", xlab = "Number of Clusters",
     main = "CRD", locate = TRUE)

plot(crt, ypar = "power", xpar = "n2",
     ylim = c(0, 1), xlim = c(10, 200), 
     ylab = "Power (ES = .20)", xlab = "Number of Clusters",
     main = "CRT", locate = TRUE)
plot(crd, ypar = "power", xpar = "n2",
     ylim = c(0, 1), xlim = c(10, 200), 
     ylab = "Power (ES = .20)", xlab = "Number of Clusters",
     main = "CRD", locate = TRUE)


