\docType{package}
\name{cosa-package}
\alias{cosa}
\alias{cosa-package}
\title{Constrained Optimal Sample Allocation}

\description{
Constrained Optimal Sample Allocation (COSA) functions are designed to optimize sample sizes at two or more levels subject to budget, statistical power, or effect size constraints. Specifying \code{rhots = 0} produces result equivalent to random assignment designs. Therefore, COSA functions also allows optimization of proportion of treatment allocation (\code{p}) under unequal marginal costs. Different starting values and algorithms may produce different results especially when sample sizes at two or more levels and \code{p} are optimized. Comparing several algorithms may faciliate decisions regarding sample sizes and \code{p}.

Designs available in \pkg{cosa} package:
\tabular{rccc}{
    \emph{Design}  \tab \emph{Total Levels} \tab \emph{Treatment Level} \cr
    \code{ird1r1}  \tab 1 \tab 1  \cr
    \code{bird2r1} \tab 2 \tab 1  \cr
    \code{bird3r1} \tab 3 \tab 1  \cr
    \code{bird4r1} \tab 4 \tab 1  \cr
    \code{crd2r2}  \tab 2 \tab 2  \cr
    \code{bcrd3r2} \tab 3 \tab 2  \cr
    \code{bcrd4r2} \tab 4 \tab 2  \cr
    \code{crd3r3}  \tab 3 \tab 3  \cr
    \code{bcrd4r3} \tab 4 \tab 3  \cr
    \code{crd4r4}  \tab 4 \tab 4
}
  \code{ird}: individual-level regression discontinuity.
  \code{bird}: block/multisite individual-level regression discontinuity.
  \code{crd}: cluster-level regression discontinuity.
  \code{bird}: block/multisite cluster-level regression discontinuity.

  Design parameters follows a sequential order. Numbers at the end of a sequential parameter refers to the corresponding level. For example \code{rho2} is the proportion of variance in the outcome between level 2 units, \code{rho3} is the proportion of variance in the outcome between level 3 units. Similiarly, \code{r21} is the proportion of the variance in the outcome explained by level 1 covariates, \code{r22} is the proportion of the variance in the outcome explained by level 2 covariates and so on. Similiar naming conventions applies to other design parameters and can be found in the documentation.

}

\author{
  Metin Bulus \email{bulusmetin@gmail.com}
}

\keyword{cosa-package}
