\name{example1}
\alias{example1}
\docType{data}
\title{
Bivariate twin data example from Classic Mx Manual
}

\description{
Data set used in some of OpenMx's examples.
}

\usage{data("example1")}
\format{
  A data frame with 400 observations on the following variables.
  \describe{
    \item{\code{IDNum}}{Twin pair ID}
    \item{\code{Zygosity}}{Zygosity of the twin pair}
    \item{\code{X1}}{X variable for twin 1}
    \item{\code{Y1}}{Y variable for twin 1}
    \item{\code{X2}}{X variable for twin 2}
    \item{\code{Y2}}{Y variable for twin 2}
  }
}

\details{
Same as \link{example2} but in wide format instead of tall.
}

\source{
Classic Mx Manual.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{
data(example1)
plot(X2 ~ X1, data=example1)
}
\keyword{datasets}
