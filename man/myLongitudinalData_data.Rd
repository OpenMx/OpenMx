\name{myLongitudinalData}
\alias{myLongitudinalData}
\docType{data}
\title{
Data for a linear latent growth curve model
}

\description{
Data set used in some of OpenMx's examples.
}

\usage{data("myLongitudinalData")}
\format{
  A data frame with 500 observations on the following variables.
  \describe{
    \item{\code{x1}}{x variable and time 1}
    \item{\code{x2}}{x variable and time 2}
    \item{\code{x3}}{x variable and time 3}
    \item{\code{x4}}{x variable and time 4}
    \item{\code{x5}}{x variable and time 5}
  }
}

\details{
Linear growth model with mean intercept around 10, and slope of about 1.5.
}

\source{
Simulated.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{
data(myLongitudinalData)

matplot(t(myLongitudinalData), type='l', lty=1)
}
\keyword{datasets}
