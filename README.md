# Constrained Optimal Sample Allocation

**cosa** implements generalized constrained optimal sample allocation (COSA) framework for multilevel regression discontinuity studies and multilevel randomized trials with continuous outcomes. The package includes functions to calculate minimum detectable effect size and statistical power. COSA functions are designed to optimize sample sizes at two or more levels subject to other constraints. 

To install the package:
```{r}
install.packages("cosa")
```

Example COSA for cluster-level regression discontinuity (three-level design, treatment at level 3):
```{r}
library(cosa)
# cost constrained - optimize n1 and n3
cosa.crd3r3(constrain = "cost", cost = 50000,
            cn1 = 5, cn2 = 10, cn3 = c(50, 20),
            es = .20, rho2 = .20, rho3 = .10,
            r21 = .20, r22 = .30, r23 = .40,
            g3 = 1, p = .40, n2 = 3)
            
# Solution converged with LBFGS 
#
# Rounded solution: 
# ----------------------------------------------- 
#   n1 [n2]  n3   [p] [cost]    es 95%lcl 95%ucl power
# 17.654    3 153 0.399  49995 0.259  0.078  0.441  0.58
# ----------------------------------------------- 
# Per unit marginal costs: 
# Level 1 treatment: 5 
# Level 1 control: 5 
# Level 2 treatment: 10 
# Level 2 control: 10 
# Level 3 treatment: 50 
# Level 3 control: 20 
# ----------------------------------------------- 
# MDES = 0.259 (with power = 80) 
# power = 0.58 (for ES = 0.2) 
# ----------------------------------------------- 
# []: point constrained (fixed) 
# <<: bound constrained 
```

Specifying `rhots = 0` produces result equivalent to random assignment design where also `p` can be optimized.

```{r}
# cost constrained - optimize n1 and p
cosa.crd3r3(rhots = 0, constrain = "cost", cost = 50000, 
            cn1 = 5, cn2 = 10, cn3 = c(50, 20),
            es = .20, rho2 = .20, rho3 = .10,
            r21 = .20, r22 = .30, r23 = .40,
            g3 = 1, n1 = 25, n2 = 3)
# Results are equivalent to random assignment design 
# Solution converged with LBFGS 
# 
# Rounded solution: 
# ----------------------------------------------- 
# [n1] [n2]  n3     p [cost]   es 95%lcl 95%ucl power
# 25    3 114 0.456  50010 0.18  0.054  0.305 0.877
# ----------------------------------------------- 
# Per unit marginal costs: 
# Level 1 treatment: 5 
# Level 1 control: 5 
# Level 2 treatment: 10 
# Level 2 control: 10 
# Level 3 treatment: 50 
# Level 3 control: 20 
# ----------------------------------------------- 
# MDES = 0.18 (with power = 80) 
# power = 0.877 (for ES = 0.2) 
# ----------------------------------------------- 
# []: point constrained (fixed) 
# <<: bound constrained 
```

<sub> Copyright &copy; 2018 Metin Bulus </sub>
