\name{crd2r2}
\alias{crd2r2}
\alias{cosa.crd2r2}
\alias{power.crd2r2}
\alias{mdes.crd2r2}
\title{Cluster-level Regression Discontinuity (Two-level Design, Discontinuity at Level 2)}

\description{
  Use \code{mdes.crd2r2()} to calculate minimum detectable effect size, \code{power.crd2r2()} to calculate statistical power, and \code{cosa.crd2r2()} for constrained optimal sample allocation. COSA can be found in the following forms; (i) under budgetary constraints given marginal costs per unit, (ii) under power or MDES contraints given marginal costs per unit, and (iii) under sample size contraints for one or more levels along with (i) or (ii).
}

\usage{
cosa.crd2r2(cn1 = 0, cn2 = 0, cost = NULL,
            n1 = NULL, n2 = NULL, p = NULL,
            n0 = c(10, 100), p0 = .50, constrain = "power",
            local.solver = c("LBFGS", "SLSQP", "MMA", "COBYLA"),
            rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
            power = .80, es = .25, alpha = .05, two.tailed = TRUE,
            rho2, g2 = 0, r21 = 0, r22 = 0)

mdes.crd2r2(power = .80, alpha = .05, two.tailed = TRUE,
            rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
            rho2, r21 = 0, r22 = 0, g2 = 0,
            p = .50, n1, n2)

power.crd2r2(es = .25, alpha = .05, two.tailed = TRUE,
             rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
             rho2, r21 = 0, r22 = 0, g2 = 0,
             p = .50, n1, n2)
}

\arguments{
  \item{cn1}{marginal cost per level 1 unit in treatment and control conditions.}
  \item{cn2}{marginal cost per level 2 unit in treatment and control conditions.}
  \item{cost}{total cost or budget.}
  \item{n1}{harmonic mean (or average) of level 1 units across level 2 units.}
  \item{n2}{number of level 2 units.}
  \item{p}{proportion of level 2 units in treatment condition.}
  \item{n0}{starting values for \code{n1, n2}.}
  \item{p0}{starting value for \code{p}}
  \item{power}{statistical power (1 - \eqn{\beta}).}
  \item{es}{effect size (Cohen's d).}
  \item{alpha}{probability of type I error (\eqn{\alpha}).}
  \item{two.tailed}{logical; \code{TRUE} for two-tailed hypothesis testing.}
  \item{constrain}{character; \code{"cost"}, \code{"power"}, or \code{"mdes"}.}
  \item{local.solver}{subset of \code{c("LBFGS", "SLSQP", "MMA", "COBYLA")}}
   \item{rhots}{correlation between the treatment and the scoring variable. Specify \code{rhots = 0} to obtain results equivalent to random assignment designs.}
  \item{k1}{left truncation point (in standard deviation units from full normal distribution mean), ignored when \code{rhots} is specified.}
  \item{k2}{right truncation point (in standard deviation units from full normal distribution mean), ignored when \code{rhots} is specified.}
  \item{dists}{character; distribution of the scoring variable, \code{"normal"} or \code{"uniform"}}
  \item{rho2}{proportion of variance in the outcome between level 2 units.}
  \item{g2}{number of covariates at level 3.}
  \item{r21}{proportion of level 1 variance in the outcome explained by level 1 covariates.}
  \item{r22}{proportion of level 2 variance in the outcome explained by level 2 covariates.}
}

\value{
  \item{parms}{list of parameters used in the function.}
  \item{exact.cosa}{exact solution, total cost, [minimum detectable] effect size and (1 - \eqn{\alpha})\% confidence limits, and statistical power.}
  \item{round.cosa}{rounded solution, total cost, [minimum detectable] effect size and (1 - \eqn{\alpha})\% confidence limits, and statistical power.}
  \item{mdes}{minimum detectable effect size and (1 - \eqn{\alpha})\% confidence limits.}
  \item{power}{statistical power (1 - \eqn{\beta})}
}

\examples{
# cost constrained
cosa.crd2r2(constrain = "cost", cost = 50000,
            cn1 = 5, cn2 = c(10, 8),
            es = .20, rho2 = .20, r21 = .20, r22 = .30,
            g2 = 1, p = .50)

# MDES
mdes.crd2r2(power = .80, rho2 = .20, r21 = .20, r22 = .30,
            g2 = 1, p = .50, n1 = 16.04, n2 = 570)

# power
power.crd2r2(es = .20, rho2 = .20, r21 = .20, r22 = .30,
             g2 = 1, p = .50, n1 = 16.04, n2 = 570)
}

\keyword{cluster-level regression discontinuity}