\name{omxCheckIdentical}
\alias{omxCheckIdentical}

\title{Exact Equality Testing Function}

\description{
   This function tests whether two objects are equal.
}

\usage{
omxCheckIdentical(a, b)
}

\arguments{
   \item{a}{the first value to compare.}
   \item{b}{the second value to compare.}
}

\details{
Performs the \sQuote{identical} comparison on the two arguments. If the two arguments are not equal, then an error will be thrown. If \sQuote{a} and \sQuote{b} are equal to each other, by default the function will print a statement informing the user the test has passed.  To turn off these print statements use \code{options("mxPrintUnitTests" = FALSE)}.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\seealso{
\code{\link{omxCheckCloseEnough}}, \code{\link{omxCheckWithinPercentError}}, \code{\link{omxCheckSetEquals}}, \code{\link{omxCheckTrue}}, \code{\link{omxCheckEquals}}
}

\examples{

omxCheckIdentical(c(1, 2, 3), c(1, 2, 3))

omxCheckIdentical(FALSE, FALSE)

# Throws an error
try(omxCheckIdentical(c(1, 2, 3), c(2, 1, 3)))

}
