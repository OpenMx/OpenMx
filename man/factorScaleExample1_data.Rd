\name{factorScaleExample1}
\alias{factorScaleExample1}
\docType{data}
\title{
Example Factor Analysis Data for Scaling the Model
}

\description{
Data set used in some of OpenMx's examples.
}

\usage{data("factorScaleExample1")}
\format{
  A data frame with 200 observations on the following variables.
  \describe{
    \item{\code{X1}}{}
    \item{\code{X2}}{}
    \item{\code{X3}}{}
    \item{\code{X4}}{}
    \item{\code{X5}}{}
    \item{\code{X6}}{}
    \item{\code{X7}}{}
    \item{\code{X8}}{}
    \item{\code{X9}}{}
    \item{\code{X10}}{}
    \item{\code{X11}}{}
    \item{\code{X12}}{}
  }
}

\details{
This appears to be a three factor model with factor 1 loading on X1-X4, factor 2 on X5-X8, and factor 3 on X9-X12.
}

\source{
Simulated
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{
data(factorScaleExample1)
round(cor(factorScaleExample1), 2)
}
\keyword{datasets}
