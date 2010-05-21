\name{omxCheckCloseEnough}
\alias{omxCheckCloseEnough}

\title{Approximate Equality Testing Function}

\description{
   This function tests whether two numeric vectors or matrixes are approximately equal to one another, within a specified threshold.
}

\usage{
omxCheckCloseEnough(a, b, epsilon = 10^(-15))
}

\arguments{
   \item{a}{a numeric vector or matrix.}
   \item{b}{a numeric vector or matrix.}
   \item{epsilon}{a non-negative tolerance threshold.}
}

\details{
Arguments \sQuote{a} and \sQuote{b} must be of the same type, ie. they must be either vectors of equal dimension or matrices of equal dimension. The two arguments are compared element-wise for approximate equality.  If the absolute value of the difference of any two values is greater than the threshold, then an error will be thrown. If \sQuote{a} and \sQuote{b} are approximately equal to each other, by default the function will print a statement informing the user the test has passed.  To turn off these print statements use \code{options("mxPrintUnitTests" = FALSE)}.
}

\seealso{
\code{\link{omxCheckWithinPercentError}}, \code{\link{omxCheckIdentical}}, \code{\link{omxCheckSetEquals}}, \code{\link{omxCheckTrue}}, \code{\link{omxCheckEquals}}
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{

omxCheckCloseEnough(c(1, 2, 3), c(1.1, 1.9 ,3.0), epsilon = 0.5)

omxCheckCloseEnough(matrix(3, 3, 3), matrix(4, 3, 3), epsilon = 2)

# Throws an error
try(omxCheckCloseEnough(c(1, 2, 3), c(1.1, 1.9 ,3.0), epsilon = 0.01))

}
