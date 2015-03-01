\name{myAutoregressiveData}
\alias{myAutoregressiveData}
\docType{data}
\title{
Example data with autoregressively related columns
}

\description{
Data set used in some of OpenMx's examples.
}

\usage{data("myAutoregressiveData")}
\format{
  A data frame with 100 observations on the following variables.
  \describe{
    \item{\code{x1}}{x variable and time 1}
    \item{\code{x2}}{x variable and time 2}
    \item{\code{x3}}{x variable and time 3}
    \item{\code{x4}}{x variable and time 4}
    \item{\code{x5}}{x variable and time 5}
  }
}

\details{
The rows are independently and identically distributed, but the columns are and auto-correlation structure.
}

\source{
Simulated.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{
data(myAutoregressiveData)
round(cor(myAutoregressiveData), 2)
# note the sub-diagonal correlations (lag 1)
#  x1-x2, x2-x3, x3-x4, x4-x5
# and the second sub-diagonal correlations (lag 2)
#  x1-x3, x2-x4, x3-x5

}
\keyword{datasets}
