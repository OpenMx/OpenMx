\name{demoOneFactor}
\alias{demoOneFactor}
\docType{data}
\title{
Demonstration data for a one factor model
}

\description{
Data set used in some of OpenMx's examples.
}

\usage{data("demoOneFactor")}
\format{
  A data frame with 500 observations on the following 5 numeric variables.
  \describe{
    \item{\code{x1}}{}
    \item{\code{x2}}{}
    \item{\code{x3}}{}
    \item{\code{x4}}{}
    \item{\code{x5}}{}
  }
}

\details{
Variables x1-x5 are typically used as indicators of the latent trait.
}

\source{
Simulated.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{
data(demoOneFactor)
cov(demoOneFactor)
cor(demoOneFactor)
}
\keyword{datasets}
