\name{jointdata}
\alias{jointdata}
\docType{data}
\title{
Joint Ordinal and continuous variables to be modeled together
}

\description{
Data set used in some of OpenMx's examples.
}

\usage{data("jointdata")}
\format{
  A data frame with 250 observations on the following variables.
  \describe{
    \item{\code{z1}}{Continuous variable}
    \item{\code{z2}}{Ordinal variable with 2 levels (0, 1)}
    \item{\code{z3}}{Continuous variable}
    \item{\code{z4}}{Ordinal variable with 4 levels (0, 1, 2, 3)}
    \item{\code{z5}}{Ordinal variable with 3 levels (0, 1, 3)}
  }
}

\details{

}

\source{
Simulated.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{
data(jointdata)
par(mfrow=c(2, 3))
h <- lapply(jointdata, hist)
par(mfrow=c(1, 1))

plot(z2 ~ z1, jointdata)

}
\keyword{datasets}
