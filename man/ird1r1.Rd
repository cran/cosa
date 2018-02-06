\name{ird1r1}
\alias{ird1r1}
\alias{mdes.ird1r1}
\alias{power.ird1r1}
\title{Simple Individual-level Regression Discontinuity}

\description{
 Use \code{mdes.ird1r1()} to calculate minimum detectable effect size, and \code{power.ird1r1()} to calculate statistical power.
}

\usage{
mdes.ird1r1(power = .80, alpha = .05, two.tailed = TRUE,
            rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
            r2 = 0, g = 0, p = .50, n)

power.ird1r1(es = .25, alpha = .05, two.tailed = TRUE,
             rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
             r2 = 0, g = 0, p = .50, n)
}

\arguments{
  \item{n}{sample size.}
  \item{p}{proportion of units in treatment condition.}
  \item{power}{statistical power (1 - \eqn{\beta}).}
  \item{es}{effect size (Cohen's d).}
  \item{alpha}{probability of type I error (\eqn{\alpha}).}
  \item{two.tailed}{logical; \code{TRUE} for two-tailed hypothesis testing.}
   \item{rhots}{correlation between the treatment and the scoring variable. Specify \code{rhots = 0} to obtain results equivalent to random assignment designs.}
  \item{k1}{left truncation point (in standard deviation units from full normal distribution mean), ignored when \code{rhots} is specified.}
  \item{k2}{right truncation point (in standard deviation units from full normal distribution mean), ignored when \code{rhots} is specified.}
  \item{dists}{character; distribution of the scoring variable, \code{"normal"} or \code{"uniform"}}
  \item{g}{number of covariates.}
  \item{r2}{proportion of variance in the outcome explained by covariates.}
}

\value{
  \item{parms}{list of parameters used in the function.}
  \item{mdes}{minimum detectable effect size and (1 - \eqn{\alpha})\% confidence limits.}
  \item{power}{statistical power (1 - \eqn{\beta})}
}

\examples{
# MDES
mdes.ird1r1(power = .80, r2 = .20, g = 1, p = .50, n = 500)

# power
power.ird1r1(es = .373, r2 = .20, g = 1, p = .50, n = 500)
}

\keyword{individual-level regression discontinuity}
