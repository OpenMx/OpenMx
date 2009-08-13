\name{omxCheckEquals}
\alias{omxCheckEquals}

\title{Exact Equality Testing Function}

\description{
   This function tests whether two objects are equal.
}

\usage{
omxCheckEquals(a, b)
}

\arguments{
   \item{a}{The first value to compare.}
   \item{b}{The second value to compare.}
}

\details{
Performs the \sQuote{==} operator on the two arguments. If the two arguments are not equal, then an error will be thrown. If \sQuote{a} and \sQuote{b} are equal to each other, by default the function will print a statement informing the user the test has passed.  To turn off these print statements use \code{options("mxPrintUnitTests" = FALSE)}.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\seealso{\itemize{
\item{\code{\link{omxCheckCloseEnough}}}{for approximate equality testing within an epsilon.}
\item{\code{\link{omxCheckWithinPercentError}}}{for approximate equality testing within a percentage.}
\item{\code{\link{omxCheckSetEquals}}}{for set equality testing.}
\item{\code{\link{omxCheckTrue}}}{for boolean equality testing.}
}}


\examples{

omxCheckEquals(c(1, 2, 3), c(1, 2, 3))

omxCheckEquals(FALSE, FALSE)

# Throws an error
try(omxCheckEquals(c(1, 2, 3), c(2, 1, 3)))

}
