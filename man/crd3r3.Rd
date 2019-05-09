\name{crd3r3}
\alias{crd3r3}
\alias{bcrd4f3}
\alias{cosa.crd3r3}
\alias{cosa.bcrd4f3}
\alias{power.crd3r3}
\alias{power.bcrd4f3}
\alias{mdes.crd3r3}
\alias{mdes.bcrd4f3}
\title{Cluster-level Regression Discontinuity (Three-level Design, Discontinuity at Level 3, w/ or w/o Strata or Fixed Blocks)}

\description{
  Use \code{mdes.crd3r3()} to calculate minimum detectable effect size, \code{power.crd3r3()} to calculate statistical power, and \code{cosa.crd3r3()} for constrained optimal sample allocation.

  If higher level strata or fixed blocks exist, use \code{mdes.bcrd4f3()} to calculate minimum detectable effect size, \code{power.bcrd4f3()} to calculate statistical power, and \code{cosa.bcrd4f3()} for constrained optimal sample allocation.
}

\usage{
cosa.crd3r3(cn1 = 0, cn2 = 0, cn3 = 0, cost = NULL,
            n1 = NULL, n2 = NULL, n3 = NULL, p = NULL,
            n0 = c(10, 3, 100 + g3), p0 = .499, constrain = "power",
            round = TRUE, max.power = FALSE,
            local.solver = c("LBFGS", "SLSQP"),
            rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
            power = .80, es = .25, alpha = .05, two.tailed = TRUE,
            rho2, rho3, g3 = 0, r21 = 0, r22 = 0, r23 = 0)

cosa.bcrd4f3(cn1 = 0, cn2 = 0, cn3 = 0, cn4 = 0, cost = NULL,
            n1 = NULL, n2 = NULL, n3 = NULL, n4 = NULL, p = NULL,
            n0 = c(10, 3, 100 + g3, 5), p0 = .499, constrain = "power",
            round = TRUE, max.power = FALSE,
            local.solver = c("LBFGS", "SLSQP"),
            rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
            power = .80, es = .25, alpha = .05, two.tailed = TRUE,
            rho2, rho3, g3 = 0, r21 = 0, r22 = 0, r23 = 0)

mdes.crd3r3(power = .80, alpha = .05, two.tailed = TRUE,
            rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
            rho2, rho3, r21 = 0, r22 = 0, r23 = 0, g3 = 0,
            p = .50, n1, n2, n3)

mdes.bcrd4f3(power = .80, alpha = .05, two.tailed = TRUE,
            rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
            rho2, rho3, r21 = 0, r22 = 0, r23 = 0, g3 = 0,
            p = .50, n1, n2, n3, n4)

power.crd3r3(es = .25, alpha = .05, two.tailed = TRUE,
             rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
             rho2, rho3, r21 = 0, r22 = 0, r23 = 0, g3 = 0,
             p = .50, n1, n2, n3)

power.bcrd4f3(es = .25, alpha = .05, two.tailed = TRUE,
             rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
             rho2, rho3, r21 = 0, r22 = 0, r23 = 0, g3 = 0,
             p = .50, n1, n2, n3, n4)
}

\arguments{
  \item{cn1}{marginal cost per level 1 unit in treatment and control conditions, e.g. \code{c(10, 5)}.}
  \item{cn2}{marginal cost per level 2 unit in treatment and control conditions, e.g. \code{c(50, 30)}.}
  \item{cn3}{marginal cost per level 3 unit in treatment and control conditions, e.g. \code{c(80, 50)}.}
  \item{cn4}{marginal cost per stratum or fixed block.}
  \item{cost}{total cost or budget.}
  \item{n1}{average number of level 1 units per level 2 unit.}
  \item{n2}{average number of level 2 units per level 3 unit.}
  \item{n3}{number of level 3 units(per stratum or block, if exists).}
  \item{n4}{number of stratum or fixed blocks.}
  \item{p}{proportion of level 3 units in treatment condition.}
  \item{n0}{vector of starting values for \code{n1, n2, n3} or \code{n1, n2, n3, n4} (positional). Starting values are replaced with averages when sample sizes are constrained by bounds.}
  \item{p0}{starting value for \code{p} when \code{rhots = 0} and \code{p = NULL}. Starting value is replaced with average when \code{p} is constrained by bounds.}
  \item{power}{statistical power (1 - \eqn{\beta}).}
  \item{es}{effect size (Cohen's d).}
  \item{alpha}{probability of type I error (\eqn{\alpha}).}
  \item{two.tailed}{logical; \code{TRUE} for two-tailed hypothesis testing.}
  \item{constrain}{character; \code{"cost"}, \code{"power"}, or \code{"mdes"}.}
  \item{round}{logical; \code{TRUE} for rounded COSA solution.}
  \item{max.power}{logical; \code{TRUE} for maximizing power instead of minimizing variance.}
  \item{local.solver}{subset of \code{c("LBFGS", "SLSQP")}}
  \item{rhots}{correlation between the treatment and the scoring variable. Specify \code{rhots = 0} to obtain results equivalent to random assignment designs.}
  \item{k1}{left truncation point (in standard deviation units from full normal distribution mean), ignored when \code{rhots} is not \code{NULL} or \code{dists = "uniform"}.}
  \item{k2}{right truncation point (in standard deviation units from full normal distribution mean), ignored when \code{rhots} is not \code{NULL} or \code{dists = "uniform"}.}
  \item{dists}{distribution of the scoring variable; \code{"normal"} or \code{"uniform"}. By default, \code{dists = "normal"} specification implies a truncated normal distribution with \code{k1 = -6} and \code{k2 = 6}.}
  \item{rho2}{proportion of variance in the outcome between level 2 units (unconditional ICC2).}
  \item{rho3}{proportion of variance in the outcome between level 3 units (unconditional ICC3).}
  \item{g3}{number of covariates at level 3.}
  \item{r21}{proportion of level 1 variance in the outcome explained by level 1 covariates.}
  \item{r22}{proportion of level 2 variance in the outcome explained by level 2 covariates.}
  \item{r23}{proportion of level 3 variance in the outcome explained by level 3 covariates.}
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
cosa.crd3r3(constrain = "cost", cost = 100000,
            cn1 = 5, cn2 = c(10, 5), cn3 = c(30, 10),
            es = .20, rho2 = .20, rho3 = .10,
            r21 = .20, r22 = .30, r23 = .40,
            g3 = 1, p = .50, n1 = 25, n2 = NULL, n3 = NULL)

# minimum detectable effect size
mdes.crd3r3(power = .80, rho2 = .20, rho3 = .10,
            r21 = .20, r22 = .30, r23 = .40,
            g3 = 1, p = .50, n1 = 19.39, n2 = 3, n3 = 300)

# statistical power
power.crd3r3(es = .20, rho2 = .20, rho3 = .10,
             r21 = .20, r22 = .30, r23 = .40,
             g3 = 1, p = .50, n1 = 19.39, n2 = 3, n3 = 300)
}
