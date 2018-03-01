\name{ird1r1}
\alias{ird1r1}
\alias{bird2f1}
\alias{cosa.ird1r1}
\alias{cosa.bird2f1}
\alias{mdes.ird1r1}
\alias{mdes.bird2f1}
\alias{power.ird1r1}
\alias{power.bird2f1}
\title{Simple Individual-level Regression Discontinuity (w/ or w/o Strata or Fixed Blocks)}

\description{
 Use \code{mdes.ird1r1()} to calculate minimum detectable effect size, \code{power.ird1r1()} to calculate statistical power, and \code{cosa.ird1r1()} for constrained optimal sample allocation.

 If higher level strata or fixed blocks exist, use \code{mdes.bird2f1()} to calculate minimum detectable effect size, \code{power.bird2f1()} to calculate statistical power, and \code{cosa.bird2f1()} for constrained optimal sample allocation.
}

\usage{
cosa.ird1r1(cn1 = 0,  cost = NULL, n1 = NULL, p = NULL,
            n0 = c(400 + g1), p0 = .50, constrain = "power",
            round = TRUE, local.solver = c("LBFGS", "SLSQP", "MMA", "COBYLA"),
            rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
            power = .80, es = .25, alpha = .05, two.tailed = TRUE,
            g1 = 0, r21 = 0)

cosa.bird2f1(cn1 = 0,  cn2 = 0, cost = NULL, n1 = NULL, n2 = NULL, p = NULL,
             n0 = c(400 + g1, 5), p0 = .50, constrain = "power",
             round = TRUE, local.solver = c("LBFGS", "SLSQP", "MMA", "COBYLA"),
             rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
             power = .80, es = .25, alpha = .05, two.tailed = TRUE,
             g1 = 0, r21 = 0)

mdes.ird1r1(power = .80, alpha = .05, two.tailed = TRUE,
            rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
            r21 = 0, g1 = 0, p = .50, n1)

mdes.bird2f1(power = .80, alpha = .05, two.tailed = TRUE,
             rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
             r21 = 0, g1 = 0, p = .50, n1, n2)

power.ird1r1(es = .25, alpha = .05, two.tailed = TRUE,
             rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
             r21 = 0, g1 = 0, p = .50, n1)

power.bird2f1(es = .25, alpha = .05, two.tailed = TRUE,
              rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
              r21 = 0, g1 = 0, p = .50, n1, n2)
}

\arguments{
  \item{cn1}{marginal cost per unit in treatment and control conditions.}
  \item{cn2}{marginal cost per stratum or fixed block.}
  \item{cost}{total cost or budget.}
  \item{n1}{sample size (per stratum or block, if exists).}
  \item{n2}{number of stratum or fixed blocks.}
  \item{p}{proportion of units in treatment condition.}
  \item{n0}{starting value for \code{n1} or \code{n1, n2}. Starting value is replaced with average when sample size is constrained by bounds.}
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
  \item{dists}{character; distribution of the scoring variable, \code{"normal"} or \code{"uniform"}. By default, \code{dists = "normal"} specification implies a truncated normal distribution with \code{k1 = -6} and \code{k2 = 6}.}
  \item{g1}{number of covariates.}
  \item{r21}{proportion of variance in the outcome explained by covariates.}
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
# cost constrained - optimize p and n1
cosa.ird1r1(rhots = 0,  constrain = "cost", cost = 1000,
            cn1 = c(10,5), p = NULL, n1 = NULL)

# minimum detectable effect size
mdes.ird1r1(power = .80, r2 = .20, g = 1, p = .50, n = 500)

# statistical power
power.ird1r1(es = .373, r2 = .20, g = 1, p = .50, n = 500)
}
