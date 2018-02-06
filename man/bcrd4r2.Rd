\name{bcrd4r2}
\alias{bcrd4r2}
\alias{cosa.bcrd4r2}
\alias{mdes.bcrd4r2}
\alias{power.bcrd4r2}
\title{Block/Multisite Cluster-level Regression Discontinuity (Four-level Design, Discontinuity at Level 2)}

\description{
 Use \code{mdes.bcrd4r2()} to calculate minimum detectable effect size, \code{power.bcrd4r2()} to calculate statistical power, and use \code{cosa.bcrd4r2()} for constrained optimal sample allocation. COSA can be found in the following forms; (i) under budgetary constraints given marginal costs per unit, (ii) under power or MDES contraints given marginal costs per unit, and (iii) under sample size contraints for one or more levels along with (i) or (ii).
}

\usage{
cosa.bcrd4r2(cn1 = 0, cn2 = 0, cn3 = 0, cn4 = 0, cost = NULL,
             n1 = NULL, n2 = NULL, n3 = NULL, n4 = NULL, p = NULL,
             n0 = c(10, 3, 100, 10), p0 = .50, constrain = "power",
             local.solver = c("LBFGS", "SLSQP", "MMA", "COBYLA"),
             rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
             power = .80, es = .25, alpha = .05, two.tailed = TRUE,
             rho2, rho3, rho4, omega3, omega4,
             g4 = 0, r21 = 0, r22 = 0, r2t3 = 0, r2t4 = 0)

mdes.bcrd4r2(power = .80, alpha = .05, two.tailed = TRUE,
             rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
             rho2, rho3, rho4, omega3, omega4,
             r21 = 0, r22 = 0, r2t3 = 0, r2t4 = 0, g4 = 0,
             p = .50, n1, n2, n3, n4)

power.bcrd4r2(es = .25, alpha = .05, two.tailed = TRUE,
              rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
              rho2, rho3, rho4, omega3, omega4,
              r21 = 0, r22 = 0, r2t3 = 0, r2t4 = 0, g4 = 0,
              p = .50, n1, n2, n3, n4)
}

\arguments{
  \item{cn1}{marginal cost per level 1 unit in treatment and control conditions.}
  \item{cn2}{marginal cost per level 2 unit in treatment and control conditions.}
  \item{cn3}{marginal cost per level 3 unit.}
  \item{cn4}{marginal cost per level 4 unit.}
  \item{cost}{total cost or budget.}
  \item{n1}{harmonic mean (or average) of level 1 units across level 2 units.}
  \item{n2}{harmonic mean (or average) of level 2 units across level 3 units.}
  \item{n3}{harmonic mean (or average) of level 3 units across level 4 units.}
  \item{n4}{number of level 4 units.}
  \item{p}{proportion of level 2 units in treatment condition.}
  \item{n0}{starting values for \code{n1, n2, n3, n4}.}
  \item{p0}{starting value for \code{p}, if \code{rhots = 0} and \code{p = NULL}}
  \item{power}{statistical power (1 - \eqn{\beta}).}
  \item{es}{effect size (Cohen's d).}
  \item{alpha}{probability of type I error (\eqn{\alpha}).}
  \item{two.tailed}{logical; \code{TRUE} for two-tailed hypothesis testing.}
  \item{constrain}{character; \code{"cost"}, \code{"power"}, or \code{"mdes"}.}
  \item{local.solver}{subset of \code{c("LBFGS", "SLSQP", "MMA", "COBYLA")}.}
  \item{rhots}{correlation between the treatment and the scoring variable. Specify \code{rhots = 0} to obtain results equivalent to random assignment designs.}
  \item{k1}{left truncation point (in standard deviation units from full normal distribution mean), ignored when \code{rhots} is specified.}
  \item{k2}{right truncation point (in standard deviation units from full normal distribution mean), ignored when \code{rhots} is specified.}
  \item{dists}{distribution of the scoring variable; \code{"normal"} or \code{"uniform"}.}
  \item{rho2}{proportion of variance in the outcome between level 2 units.}
  \item{rho3}{proportion of variance in the outcome between level 3 units.}
  \item{rho4}{proportion of variance in the outcome between level 4 units.}
  \item{omega3}{ratio of the treatment effect variance between level 3 units to the variance in the outcome between level 3 units.}
  \item{omega4}{ratio of the treatment effect variance between level 4 units to the variance in the outcome between level 4 units.}
  \item{g4}{number of covariates at level 4.}
  \item{r21}{proportion of level 1 variance in the outcome explained by level 1 covariates.}
  \item{r22}{proportion of level 1 variance in the outcome explained by level 2 covariates.}
  \item{r2t3}{proportion of treatment effect variance between level 3 units explained by level 3 covariates.}
  \item{r2t4}{proportion of treatment effect variance between level 4 units explained by level 4 covariates.}
}

\value{
  \item{parms}{list of parameters used in the function.}
  \item{exact.cosa}{exact solution, total cost, [minimum detectable] effect size and (1 - \eqn{\alpha})\% confidence limits, and statistical power.}
  \item{round.cosa}{rounded solution, total cost, [minimum detectable] effect size and (1 - \eqn{\alpha})\% confidence limits, and statistical power.}
  \item{mdes}{minimum detectable effect size and (1 - \eqn{\alpha})\% confidence limits.}
  \item{power}{statistical power (1 - \eqn{\beta})}
}

\examples{
  \donttest{
# power constrained
cosa.bcrd4r2(constrain = "power", power = .80,
             es = .20, rho2 = .20, rho3 = .10, rho4 = .05,
             omega3 = .10, omega4 = .20,
             r21 = .20, r22 = .20, r2t3 = .20, r2t4 = .20,
             g4 = 1, p = .50,
             n1 = c(10, 40), n2 = 3, n4 = 10)
  }

# cost constrained
cosa.bcrd4r2(constrain = "cost", cost = 100000,
             cn1 = c(5, 3), cn2 = 10, cn3 = 30, cn4 = 50,
             es = .20, rho2 = .20, rho3 = .10, rho4 = .05,
             omega3 = .10, omega4 = .20,
             r21 = .20, r22 = .20, r2t3 = .20, r2t4 = .20,
             g4 = 1, p = .50,
             n1 = 25, n2 = 3)

# MDES
mdes.bcrd4r2(power = .80, rho2 = .20, rho3 = .10, rho4 = .05,
             omega3 = .10, omega4 = .20,
             r21 = .20, r22 = .30, r2t3 = .40, r2t4 = .50,
             g4 = 1, p = .50, n1 = 10, n2 = 3, n3 = 50.92, n4 = 13)

# power
power.bcrd4r2(es = .20, rho2 = .20, rho3 = .10, rho4 = .05,
              omega3 = .10, omega4 = .20,
              r21 = .20, r22 = .30, r2t3 = .40, r2t4 = .50,
              g4 = 1, p = .50, n1 = 10, n2 = 3, n3 = 50.92, n4 = 13)
}

\keyword{cluster-level regression discontinuity}
