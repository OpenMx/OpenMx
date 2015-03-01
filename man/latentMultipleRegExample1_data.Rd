\name{latentMultipleRegExample1}
\alias{latentMultipleRegExample1}
\docType{data}
\title{
Example data for multiple regression among latent variables
}

\description{
Data set used in some of OpenMx's examples.
}

\usage{data("latentMultipleRegExample1")}
\format{
  A data frame with 200 observations on the following variables.
  \describe{
    \item{\code{X1}}{Factor 1 indicator}
    \item{\code{X2}}{Factor 1 indicator}
    \item{\code{X3}}{Factor 1 indicator}
    \item{\code{X4}}{Factor 1 indicator}
    \item{\code{X5}}{Factor 2 indicator}
    \item{\code{X6}}{Factor 2 indicator}
    \item{\code{X7}}{Factor 2 indicator}
    \item{\code{X8}}{Factor 2 indicator}
    \item{\code{X9}}{Factor 3 indicator}
    \item{\code{X10}}{Factor 3 indicator}
    \item{\code{X11}}{Factor 3 indicator}
    \item{\code{X12}}{Factor 3 indicator}
  }
}

\details{
Factor 1 strongly predicts factor 3.  Factor 2 weakly predicts factor 3.
}

\source{
Simulated.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{
data(latentMultipleRegExample1)
round(cor(latentMultipleRegExample1), 2)
}
\keyword{datasets}
