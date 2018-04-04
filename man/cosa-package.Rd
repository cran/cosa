\docType{package}
\name{cosa-package}
\alias{cosa-package}
\alias{cosa}
\alias{mdes}
\alias{power}
\title{Constrained Optimal Sample Allocation}

\description{
Constrained Optimal Sample Allocation (COSA) functions are designed to optimize sample sizes at one or more levels subject to budget, statistical power, or effect size constraints. COSA can be found in the following forms; (i) under budgetary constraints given marginal costs per unit while minimizing the sampling variance of the treatment effect (or, alternatively, while maximizing statistical power), (ii) under statistical power or effect size (ES) constraints while minimizing the total cost, and (iii) under sample size constraints for one or more levels along with (i) or (ii). Specifying \code{rhots = 0} produces result equivalent to corresponding random assignment designs, which means there is no relationship between the treatment [random] and the score variable. Therefore, COSA functions also allows optimization of proportion of treatment allocation (\code{p}) under unequal marginal costs when primary constraint is on total cost. Different starting values and algorithms may produce different results especially when marginal cost information is not provided and sample sizes at two or more levels and \code{p} are optimized. Comparing several algorithms may faciliate decisions regarding sample sizes and \code{p}.

Designs available in \pkg{cosa} package:
\tabular{rccc}{
    \emph{Design}  \tab \emph{Total Levels} \tab \emph{Treatment Level}  \tab \emph{Top Level}\cr
    \code{ird1r1}  \tab 1 \tab 1 \tab random  \cr
    \code{bird2r1} \tab 2 \tab 1 \tab random  \cr
    \code{bird2f1} \tab 2 \tab 1 \tab fixed   \cr
    \code{bird3r1} \tab 3 \tab 1 \tab random  \cr
    \code{bird4r1} \tab 4 \tab 1 \tab random  \cr
    \code{crd2r2}  \tab 2 \tab 2 \tab random  \cr
    \code{bcrd3f2} \tab 3 \tab 2 \tab fixed   \cr
    \code{bcrd3r2} \tab 3 \tab 2 \tab random  \cr
    \code{bcrd4r2} \tab 4 \tab 2 \tab random  \cr
    \code{crd3r3}  \tab 3 \tab 3 \tab random  \cr
    \code{bcrd4f3} \tab 4 \tab 3 \tab fixed   \cr
    \code{bcrd4r3} \tab 4 \tab 3 \tab random  \cr
    \code{crd4r4}  \tab 4 \tab 4 \tab random
}
  \code{ird}: individual-level regression discontinuity.
  \code{bird}: blocked individual-level regression discontinuity.
  \code{crd}: cluster-level regression discontinuity.
  \code{bcrd}: blocked cluster-level regression discontinuity.

  Design parameters follows a sequential order. Numbers at the end of a sequential parameter refers to the corresponding level. For example \code{rho2} is the proportion of variance in the outcome between level 2 units, \code{rho3} is the proportion of variance in the outcome between level 3 units. Similiarly, \code{r21} is the proportion of the variance in the outcome explained by level 1 covariates, \code{r22} is the proportion of the variance in the outcome explained by level 2 covariates and so on. Similiar naming conventions applies to other design parameters.

  Wrapper functions described below pass arguments to design-specific functions. Most of the arguments should be specified in the vector or list form length of which depends on the type of design specified in \code{design} argument. To use these functions, some knowledge about design parameters is required. For example for \code{design = "bird3r1"} the length of \code{rho} is 2, the length of \code{omega} is 2, the length of \code{r2} is 1, and the length of \code{r2t} is 2. Otherwise design-specific functions are recommended.

}

\usage{
cosa(design = "ird1r1",
     cn = as.list(rep(0, length(n))), cost = NULL,
     n = list(NULL), p = NULL, n0 = rep(5 + g, length(n)), p0 = .499,
     constrain = "power", local.solver = c("LBFGS", "SLSQP", "MMA", "COBYLA"),
     power = .80, es = .25, alpha = .05, two.tailed = TRUE,
     rhots = NULL, k = c(-6, 6), dists = "normal",
     rho = NULL, omega = NULL, r2t = NULL,
     r2 = 0,  g = 0)

mdes(design = "ird1r1", power = .80, alpha = .05, two.tailed = TRUE,
     rhots = NULL, k = c(-6, 6), dists = "normal",
     rho = NULL, omega = NULL, r2t = NULL,
     r2 = 0, n = 250, g = 0, p = .50)

power(design = "ird1r1", es = .25, alpha = .05, two.tailed = TRUE,
      rhots = NULL, k = c(-6, 6), dists = "normal",
      rho = NULL, omega = NULL, r2t = NULL,
      r2 = 0, n = 250, g = 0, p = .50)
}

\arguments{
  \item{design}{One of the \code{"ird1r1"}, \code{"bird2f1"}, \code{"bird2r1"}, \code{"bird3r1"}, \code{"bird4r1"}, \code{"crd2r2"}, \code{"bcrd3f2"}, \code{"bcrd3r2"}, \code{"bcrd4r2"}, \code{"crd3r3"}, \code{"bcrd4f3"}, \code{"bcrd4r3"}, \code{"crd4r4"}.}
  \item{cn}{list of marginal costs per unit (positional). e.g. \code{list(cn1, cn2, cn3)}.}
  \item{cost}{total cost or budget.}
  \item{n}{list of sample sizes (positional). e.g. \code{list(n1, n2, n3)}.}
  \item{p}{proportion of units in treatment condition.}
  \item{n0}{vector of starting values for \code{n} (positional). Starting values are replaced with averages when sample sizes are constrained by bounds.}
  \item{p0}{starting value for \code{p} when \code{rhots = 0} and \code{p = NULL}. Starting value is replaced with average when \code{p} is constrained by bounds.}
  \item{power}{statistical power (1 - \eqn{\beta}).}
  \item{es}{effect size (Cohen's d).}
  \item{alpha}{probability of type I error (\eqn{\alpha}).}
  \item{two.tailed}{logical; \code{TRUE} for two-tailed hypothesis testing.}
  \item{constrain}{character; \code{"cost"}, \code{"power"}, or \code{"mdes"}.}
  \item{local.solver}{subset of \code{c("LBFGS", "SLSQP", "MMA", "COBYLA")}.}
  \item{rhots}{correlation between the treatment and the scoring variable. Specify \code{rhots = 0} to obtain results equivalent to random assignment designs.}
  \item{k}{vector of left and right truncation point (in standard deviation units from full normal distribution mean), ignored when \code{rhots} is not \code{NULL} or \code{dists = "uniform"}. e.g. \code{c(k1, k2)}. }
  \item{dists}{character; distribution of the scoring variable, \code{"normal"} or \code{"uniform"}. By default, \code{dists = "normal"} specification implies a truncated normal distribution with \code{k = c(-6, 6)}.}
  \item{rho}{vector or list; proportion of variance in the outcome between level X units. e.g. \code{c(rho2, rho3)} (unconditional ICC).}
  \item{omega}{vector or list; ratio of the treatment effect variance between level X units to the variance in the outcome between level X units (positional) e.g. \code{c(omega2, omega3)}.}
  \item{g}{number of covariates at the top level.}
  \item{r2}{vector or list; proportion of  variance in the outcome explained by level X covariates (positional). e.g. \code{c(r21, r23)}}
  \item{r2t}{vector or list; proportion of treatment effect variance between level X units explained by level X covariates (positional). e.g. \code{c(r2t2, r2t3)}}
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
cosa(design = "crd3r3", rhots = 0,
     cn =list(c(5, 3), c(10,5), c(50,20)),
     constrain = "cost", cost = 15000,
     rho = c(.17, .06),
     g = 4, r2 = c(.30, .20, .20),
     n = list(15, NULL, NULL))

# minimum detectable effect size
mdes(design = "crd2r2",
     rho = .20,
     r2 = c(0, 0),
     p = .50, n = c(4, 20))

# statistical power
power(design = "bird3r1",
      rho = c(.17, .06),
      r2 = 0, r2t = c(0, 0),
      omega = c(.30, .30),
      p = .50, n = c(15, 3, 60))
}


\author{
  Metin Bulus \email{bulusmetin@gmail.com}
}
