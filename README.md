# Bound Constrained Optimal Sample Allocation

**cosa** implements (generalized) bound constrained optimal sample allocation framework for multilevel regression discontinuity designs (RDDs) and multilevel randomized controlled trials (CRTs) with continuous outcomes. COSA functions are designed to optimize proportion of treatment allocation (`p`) and sample size (generically, `n`) at one or more levels subject to budget, statistical power, or effect size constraints along with constraints on `n` and `p`. Constraints on `n` and `p` can be in the form of fixed values or bound constraints (aka box constraints).

Note: `n` and `p` should be omitted (or specified as `NULL`) for optimization. `p` can only be optimized (or bound constrained) in multilevel randomized trials when treatment and control units have differing costs. When the primary constraint is on statistical power or effect size providing marginal cost information will result in an cost-efficient sample allocation. When marginal cost information is not provided different starting values and algorithms may produce different results especially when sample sizes at two or more levels and `p` are optimized. Comparing several algorithms and starting values may faciliate decisions regarding sample sizes and `p`.

To install and load the package:
```{r}
install.packages("cosa")
library(cosa)
```

Constrained optimal sample allocation in cluster-level regression discontinuity study where discontinuity resides at level 3:

```{r}
# cost constrained - optimize n1 and n3
rdd <- cosa.crd3r3(constrain = "cost", cost = 30000,
                   cn1 = 5, cn2 = 10, cn3 = c(50, 20),
                   es = .20, rho2 = .20, rho3 = .10,
                   g3 = 1, r21 = .20, r22 = .30, r23 = .40,
                   p = .40, n1 = NULL, n2 = 3, n3 = NULL)
summary(rdd)
plot(rdd)
```

Specifying `rhots = 0` produces result equivalent to corresponding random assignment design where also `p` can be optimized. `rhots = 0` means there is no relationship between the treatment [random] and the score variable. Constrained optimal sample allocation in three-level cluster randomized trial where treatment resides at level 3:

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
