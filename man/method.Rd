\name{summary}
\alias{summary.power}
\alias{summary.mdes}
\alias{summary.cosa}
\title{Summary}

\description{
  Prints a summary of a \pkg{cosa} object.
}

\usage{
  \method{summary}{power}(object, ...)
  \method{summary}{mdes}(object, ...)
  \method{summary}{cosa}(object, ...)
}

\arguments{
  \item{object}{an object returned from one of the functions in \pkg{cosa} package.}
  \item{...}{for S3 consistency, has no effect on summary.}
}


\examples{
# cost constrained
d1 <- cosa.bcrd3r2(constrain = "cost", cost = 20000,
             cn1 = c(5, 3), cn2 = 10, cn3 = 30,
             es = .20, rho2 = .20, rho3 = .10, omega3 = .10,
             r21 = .20, r22 = .20, r2t3 = .20,
             g3 = 1, p = .50, n1 = 25)
summary(d1)
}
