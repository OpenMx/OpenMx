\name{multiData1}
\alias{multiData1}
\docType{data}
\title{
Data for multiple regression
}

\description{
Data set used in some of OpenMx's examples.
}

\usage{data("multiData1")}
\format{
  A data frame with 500 observations on the following variables.
  \describe{
    \item{\code{x1}}{}
    \item{\code{x2}}{}
    \item{\code{x3}}{}
    \item{\code{x4}}{}
    \item{\code{y}}{}
  }
}

\details{
x1-x4 are predictor variables, and y is the outcome.
}

\source{
Simulated.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{
data(multiData1)
summary(lm(y ~ ., data=multiData1))
#results can be replicated in OpenMx.

}
\keyword{datasets}
