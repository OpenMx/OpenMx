\name{omxCheckSetEquals}
\alias{omxCheckSetEquals}

\title{Set Equality Testing Function}

\description{
   This function tests whether two vectors contain the same elements.
}

\usage{
omxCheckSetEquals(a, b)
}

\arguments{
   \item{a}{the first vector to compare.}
   \item{b}{the second vector to compare.}
}

\details{
Performs the \sQuote{setequal} function on the two arguments. If the two arguments do not contain the same elements, then an error will be thrown. If \sQuote{a} and \sQuote{b} contain the same elements, by default the function will print a statement informing the user the test has passed.  To turn off these print statements use \code{options("mxPrintUnitTests" = FALSE)}.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\seealso{
\code{\link{omxCheckCloseEnough}}, \code{\link{omxCheckWithinPercentError}}, \code{\link{omxCheckIdentical}}, \code{\link{omxCheckTrue}}, \code{\link{omxCheckEquals}}
}

\examples{

omxCheckSetEquals(c(1, 1, 2, 2, 3), c(3, 2, 1))

omxCheckSetEquals(matrix(1, 1, 1), matrix(1, 3, 3))

# Throws an error
try(omxCheckSetEquals(c(1, 2, 3, 4), c(2, 1, 3)))

}
