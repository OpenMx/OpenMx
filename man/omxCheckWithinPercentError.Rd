\name{omxCheckWithinPercentError}
\alias{omxCheckWithinPercentError}

\title{Approximate Percent Equality Testing Function}

\description{
   This function tests whether two numeric vectors or matrixes are approximately equal to one another, within a specified percentage.
}

\usage{
omxCheckWithinPercentError(a, b, percent = 0.1)
}

\arguments{
   \item{a}{A numeric vector or matrix.}
   \item{b}{A numeric vector or matrix.}
   \item{percent}{A non-negative percentage.}
}

\details{
Arguments \sQuote{a} and \sQuote{b} must be of the same type, ie. they must be either vectors of equal dimension or matrices of equal dimension. The two arguments are compared element-wise for approximate equality.  If the absolute value of the difference of any two values is greater than the percentage difference of \sQuote{a}, then an error will be thrown. If \sQuote{a} and \sQuote{b} are approximately equal to each other, by default the function will print a statement informing the user the test has passed.  To turn off these print statements use \code{options("mxPrintUnitTests" = FALSE)}.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\seealso{\itemize{
\item{\code{\link{omxCheckCloseEnough}}}{for approximate equality testing within an epsilon.}
\item{\code{\link{omxCheckIdentical}}}{for exact equality testing.}
\item{\code{\link{omxCheckSetEquals}}}{for set equality testing.}
\item{\code{\link{omxCheckTrue}}}{for boolean equality testing.}
}}


\examples{

omxCheckWithinPercentError(c(1, 2, 3), c(1.1, 1.9 ,3.0), percent = 50)

omxCheckWithinPercentError(matrix(3, 3, 3), matrix(4, 3, 3), percent = 150)

# Throws an error
try(omxCheckWithinPercentError(c(1, 2, 3), c(1.1, 1.9 ,3.0), percent = 0.01))

}
