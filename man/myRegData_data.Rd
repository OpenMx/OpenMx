\name{myRegData}
\alias{myRegData}
\docType{data}
\title{
Example regression data with correlated predictors
}

\description{
Data set used in some of OpenMx's examples.
}

\usage{data("myRegData")}
\format{
  A data frame with 100 observations on the following variables.
  \describe{
    \item{\code{w}}{Predictor variable}
    \item{\code{x}}{Predictor variable}
    \item{\code{y}}{Predictor variable}
    \item{\code{z}}{Outcome varialbe}
  }
}

\details{
w, x, and y are predictors of z.  x and y are correlated.
}

\source{
Simulated.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{
data(myRegData)
summary(lm(z ~ ., data=myRegData))

}
\keyword{datasets}
