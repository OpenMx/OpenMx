\name{demoTwoFactor}
\alias{demoTwoFactor}
\docType{data}
\title{
Demonstration data for a two factor model
}

\description{
Data set used in some of OpenMx's examples.
}

\usage{data("demoTwoFactor")}
\format{
  A data frame with 500 observations on the following 10 numeric variables.
  \describe{
    \item{\code{x1}}{}
    \item{\code{x2}}{}
    \item{\code{x3}}{}
    \item{\code{x4}}{}
    \item{\code{x5}}{}
    \item{\code{y1}}{}
    \item{\code{y2}}{}
    \item{\code{y3}}{}
    \item{\code{y4}}{}
    \item{\code{y5}}{}
  }
}

\details{
Variables x1-x5 are typically used as indicators of one latent trait.  Variables y1-y5 are typically used as indicators of another latent trait.
}

\source{
Simulated.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{
data(demoTwoFactor)
cov(demoTwoFactor)
cor(demoTwoFactor)
}
\keyword{datasets}
