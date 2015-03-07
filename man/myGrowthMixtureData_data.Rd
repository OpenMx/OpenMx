\name{myGrowthMixtureData}
\alias{myGrowthMixtureData}
\docType{data}
\title{
Data for a growth mixture model
}

\description{
Data set used in some of OpenMx's examples.
}

\usage{data("myGrowthMixtureData")}
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
The same as \link{myGrowthKnownClassData}, but without the class membership variable.
}

\source{
Simulated.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{
data(myGrowthMixtureData)

matplot(t(myGrowthMixtureData), type='l', lty=1)

data(myGrowthKnownClassData)
all(myGrowthKnownClassData[,-6]==myGrowthMixtureData)
}
\keyword{datasets}
