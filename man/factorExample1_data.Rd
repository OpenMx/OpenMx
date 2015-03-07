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
    \item{\code{x1}}{}
    \item{\code{x2}}{}
    \item{\code{x3}}{}
    \item{\code{x4}}{}
    \item{\code{x5}}{}
    \item{\code{x6}}{}
    \item{\code{x7}}{}
    \item{\code{x8}}{}
    \item{\code{x9}}{}
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
