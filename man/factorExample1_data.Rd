\name{factorExample1}
\alias{factorExample1}
\docType{data}
\title{
Example Factor Analysis Data
}

\description{
Data set used in some of OpenMx's examples.
}

\usage{data("factorExample1")}
\format{
  A data frame with 500 observations on the following variables.
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
  }
}

\details{
This appears to be a three factor model, but perhaps with an odd loading structure.
}

\source{
Simulated
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{
data(factorExample1)
round(cor(factorExample1), 2)

factanal(covmat=cov(factorExample1), factors=3, rotation="promax")

}
\keyword{datasets}
