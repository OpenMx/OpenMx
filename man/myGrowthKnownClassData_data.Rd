\name{myGrowthKnownClassData}
\alias{myGrowthKnownClassData}
\docType{data}
\title{
Data for a growth mixture model with the true class membership
}

\description{
Data set used in some of OpenMx's examples.
}

\usage{data("myGrowthKnownClassData")}
\format{
  A data frame with 500 observations on the following variables.
  \describe{
    \item{\code{x1}}{x variable and time 1}
    \item{\code{x2}}{x variable and time 2}
    \item{\code{x3}}{x variable and time 3}
    \item{\code{x4}}{x variable and time 4}
    \item{\code{x5}}{x variable and time 5}
    \item{\code{c}}{Known class membership variable}
  }
}

\details{
The same as \link{myGrowthMixtureData}, but with the class membership variable.
}

\source{
Simulated.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{
data(myGrowthKnownClassData)

#plot the observed trajectories
# blue lines are class 1, green lines are class 2
colSel <-c('blue', 'green')[myGrowthKnownClassData$c]
matplot(t(myGrowthKnownClassData[,-6]), type='l', lty=1, col=colSel)

}
\keyword{datasets}
