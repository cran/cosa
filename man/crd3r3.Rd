\name{crd3r3}
\alias{crd3r3}
\alias{cosa.crd3r3}
\alias{power.crd3r3}
\alias{mdes.crd3r3}
\title{Cluster-level Regression Discontinuity (Three-level Design, Discontinuity at Level 3)}

\description{
  Use \code{mdes.crd3r3()} to calculate minimum detectable effect size, \code{power.crd3r3()} to calculate statistical power, and \code{cosa.crd3r3()} for constrained optimal sample allocation. COSA can be found in the following forms; (i) under budgetary constraints given marginal costs per unit, (ii) under power or MDES contraints given marginal costs per unit, and (iii) under sample size contraints for one or more levels along with (i) or (ii).
}

\usage{
cosa.crd3r3(cn1 = 0, cn2 = 0, cn3 = 0, cost = NULL,
            n1 = NULL, n2 = NULL, n3 = NULL, p = NULL,
            n0 = c(10, 3, 100), p0 = .50, constrain = "power",
            local.solver = c("LBFGS", "SLSQP", "MMA", "COBYLA"),
            rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
            power = .80, es = .25, alpha = .05, two.tailed = TRUE,
            rho2, rho3, g3 = 0, r21 = 0, r22 = 0, r23 = 0)

mdes.crd3r3(power = .80, alpha = .05, two.tailed = TRUE,
            rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
            rho2, rho3, r21 = 0, r22 = 0, r23 = 0, g3 = 0,
            p = .50, n1, n2, n3)

power.crd3r3(es = .25, alpha = .05, two.tailed = TRUE,
             rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
             rho2, rho3, r21 = 0, r22 = 0, r23 = 0, g3 = 0,
             p = .50, n1, n2, n3)
}

\arguments{
  \item{cn1}{marginal cost per level 1 unit in treatment and control conditions.}
  \item{cn2}{marginal cost per level 2 unit in treatment and control conditions.}
  \item{cn3}{marginal cost per level 3 unit in treatment and control conditions.}
  \item{cost}{total cost or budget.}
  \item{n1}{harmonic mean (or average) of level 1 units across level 2 units.}
  \item{n2}{harmonic mean (or average) of level 2 units across level 3 units.}
  \item{n3}{number of level 3 units across level 4 units.}
  \item{p}{proportion of level 3 units in treatment condition.}
  \item{n0}{starting values for \code{n1, n2, n3}.}
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
  \item{dists}{distribution of the scoring variable; \code{"normal"} or \code{"uniform"}.}
  \item{rho2}{proportion of variance in the outcome between level 2 units.}
  \item{rho3}{proportion of variance in the outcome between level 3 units.}
  \item{g3}{number of covariates at level 3.}
  \item{r21}{proportion of level 1 variance in the outcome explained by level 1 covariates.}
  \item{r22}{proportion of level 2 variance in the outcome explained by level 2 covariates.}
  \item{r23}{proportion of level 3 variance in the outcome explained by level 3 covariates.}
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
cosa.crd3r3(constrain = "cost", cost = 100000,
            cn1 = 5, cn2 = c(10, 5), cn3 = c(30, 10),
            es = .20, rho2 = .20, rho3 = .10,
            r21 = .20, r22 = .30, r23 = .40,
            g3 = 1, p = .50,
            n1 = 25)

# MDES
mdes.crd3r3(power = .80, rho2 = .20, rho3 = .10,
            r21 = .20, r22 = .30, r23 = .40,
            g3 = 1, p = .50, n1 = 19.39, n2 = 3, n3 = 300)

# power
power.crd3r3(es = .20, rho2 = .20, rho3 = .10,
             r21 = .20, r22 = .30, r23 = .40,
             g3 = 1, p = .50, n1 = 19.39, n2 = 3, n3 = 300)
}

\keyword{cluster-level regression discontinuity}
