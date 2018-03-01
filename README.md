# Constrained Optimal Sample Allocation

**cosa** implements generalized constrained optimal sample allocation (COSA) framework for multilevel regression discontinuity studies and multilevel randomized trials with continuous outcomes. COSA functions are designed to optimize proportion of treatment allocation (`p`) and sample sizes (`n`) at one or more levels subject to budget, statistical power, or effect size constraints along with constraints on `n` and `p`. Constraints on `n` and `p` can be in the form of fixed values or bound constraints (or box constraints).

Note: `n` and `p` should be omitted (or specified as `NULL`) for optimization. `p` can only be optimized (or bound constrained) in multilevel randomized trials when treatment and control units have differing costs and when primary constraint is on the total cost. 

To install the package:
```{r}
install.packages("cosa")
```

Example COSA for cluster-level regression discontinuity (three-level design, treatment at level 3):
```{r}
library(cosa)
# cost constrained - optimize n1 and n3
rdd <- cosa.crd3r3(constrain = "cost", cost = 30000,
                   cn1 = 5, cn2 = 10, cn3 = c(50, 20),
                   es = .20, rho2 = .20, rho3 = .10,
                   g3 = 1, r21 = .20, r22 = .30, r23 = .40,
                   p = .40, n1 = NULL, n2 = 3, n3 = NULL)
summary(rdd)
plot(rdd)
```

Specifying `rhots = 0` produces result equivalent to corresponding random assignment design where also `p` can be optimized.

```{r}
# cost constrained - optimize p, n1, n3
crt <- cosa.crd3r3(rhots = 0, constrain = "cost", cost = 30000, 
                   cn1 = 5, cn2 = 10, cn3 = c(50, 20),
                   es = .20, rho2 = .20, rho3 = .10,
                   g3 = 1, r21 = .20, r22 = .30, r23 = .40,
                   p = NULL, n1 = NULL, n2 = 3, n3 = NULL)
summary(crt)
plot(crt)
```

--o-- <br>
<footnote> Copyright &copy; 2017-2018 Metin Bulus </footnote>
