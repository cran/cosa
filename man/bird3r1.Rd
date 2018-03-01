\name{bird3r1}
\alias{bird3r1}
\alias{cosa.bird3r1}
\alias{mdes.bird3r1}
\alias{power.bird3r1}
\title{Blocked Individual-level Regression Discontinuity (Three-level Design, Discontinuity at Level 1)}

\description{
 Use \code{mdes.bird3r1()} to calculate minimum detectable effect size, \code{power.bird3r1()} to calculate statistical power, and \code{cosa.bird3r1()} for constrained optimal sample allocation.
}

\usage{
cosa.bird3r1(cn1 = 0, cn2 = 0, cn3 = 0, cost = NULL,
             n1 = NULL, n2 = NULL, n3 = NULL, p = NULL,
             n0 = c(10, 3, 100 + g3), p0 = .50, constrain = "power",
             round = TRUE, local.solver = c("LBFGS", "SLSQP", "MMA", "COBYLA"),
             rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
             power = .80, es = .25, alpha = .05, two.tailed = TRUE,
             rho2, rho3, omega2, omega3,
             g3 = 0, r21 = 0, r2t2 = 0, r2t3 = 0)

mdes.bird3r1(power = .80, alpha = .05, two.tailed = TRUE,
             rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
             rho2, rho3, omega2, omega3,
             r21 = 0, r2t2 = 0, r2t3 = 0, g3 = 0,
             p = .50, n1, n2, n3)

power.bird3r1(es = .25, alpha = .05, two.tailed = TRUE,
              rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
              rho2, rho3, omega2, omega3,
              r21 = 0, r2t2 = 0, r2t3 = 0, g3 = 0,
              p = .50, n1, n2, n3)
}

\arguments{
  \item{cn1}{marginal cost per level 1 unit in treatment and control conditions.}
  \item{cn2}{marginal cost per level 2 unit.}
  \item{cn3}{marginal cost per level 3 unit.}
  \item{cost}{total cost or budget.}
  \item{n1}{average number of level 1 units per level 2 unit.}
  \item{n2}{average number of level 2 units per level 3 unit.}
  \item{n3}{number of level 3 units.}
  \item{p}{proportion of level 1 units in treatment condition.}
  \item{n0}{vector of starting values for \code{n1, n2, n3} (positional). Starting values are replaced with averages when sample sizes are constrained by bounds.}
  \item{p0}{starting value for \code{p} when \code{rhots = 0} and \code{p = NULL}. Starting value is replaced with average when \code{p} is constrained by bounds.}
  \item{power}{statistical power (1 - \eqn{\beta}).}
  \item{es}{effect size (Cohen's d).}
  \item{alpha}{probability of type I error (\eqn{\alpha}).}
  \item{two.tailed}{logical; \code{TRUE} for two-tailed hypothesis testing.}
  \item{constrain}{character; \code{"cost"}, \code{"power"}, or \code{"mdes"}.}
  \item{round}{logical; \code{TRUE} for rounded COSA solution.}
  \item{local.solver}{subset of \code{c("LBFGS", "SLSQP", "MMA", "COBYLA")}}
  \item{rhots}{correlation between the treatment and the scoring variable. Specify \code{rhots = 0} to obtain results equivalent to random assignment designs.}
  \item{k1}{left truncation point (in standard deviation units from full normal distribution mean), ignored when \code{rhots} is not \code{NULL} or \code{dists = "uniform"}.}
  \item{k2}{right truncation point (in standard deviation units from full normal distribution mean), ignored when \code{rhots} is not \code{NULL} or \code{dists = "uniform"}.}
  \item{dists}{distribution of the scoring variable; \code{"normal"} or \code{"uniform"}. By default, \code{dists = "normal"} specification implies a truncated normal distribution with \code{k1 = -6} and \code{k2 = 6}.}
  \item{rho2}{proportion of variance in the outcome between level 2 units (unconditional ICC2).}
  \item{rho3}{proportion of variance in the outcome between level 3 units (unconditional ICC3).}
  \item{omega2}{ratio of the treatment effect variance between level 2 units to the variance in the outcome between level 2 units.}
  \item{omega3}{ratio of the treatment effect variance between level 3 units to the variance in the outcome between level 3 units.}
  \item{g3}{number of covariates at level 3.}
  \item{r21}{proportion of level 1 variance in the outcome explained by level 1 covariates.}
  \item{r2t2}{proportion of treatment effect variance between level 2 units explained by level 2 covariates.}
  \item{r2t3}{proportion of treatment effect variance between level 3 units explained by level 3 covariates.}
}

\value{
  \item{parms}{list of parameters used in the function.}
  \item{df}{degrees of freedom.}
  \item{sse}{standardized standard error.}
  \item{cosa}{constrained optimal sample allocation.}
  \item{mdes}{minimum detectable effect size and (1 - \eqn{\alpha})\% confidence limits.}
  \item{power}{statistical power (1 - \eqn{\beta})}
}

\examples{
# cost constrained - optimize n2 and n3
cosa.bird3r1(constrain = "cost", cost = 5000,
             cn1 = c(5, 3), cn2 = 10, cn3 = 30,
             es = .20, rho2 = .20, rho3 = .10, omega2 = .20, omega3 = .10,
             r21 = .20, r2t2 = .20, r2t3 = .20, g3 = 1, p = .50,
             n1 = 25, n2 = NULL, n3 = NULL)

# minimum detectable effect size
mdes.bird3r1(power = .80, rho2 = .20, rho3 = .10, omega2 = .20, omega3 = .10,
             r21 = .20, r2t2 = .30, r2t3 = .40, g3 = 1, p = .50,
             n1 = 28.06, n2 = 3, n3 = 16)

# statistical power
power.bird3r1(es = .20, rho2 = .20, rho3 = .10, omega2 = .20, omega3 = .10,
              r21 = .20, r2t2 = .30, r2t3 = .40, g3 = 1, p = .50,
              n1 = 28.06, n2 = 3, n3 = 16)

}
