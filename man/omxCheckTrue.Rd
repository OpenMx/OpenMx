\name{omxCheckTrue}
\alias{omxCheckTrue}

\title{Boolean Equality Testing Function}

\description{
   This function tests whether an object is equal to TRUE.
}

\usage{
omxCheckTrue(a)
}

\arguments{
   \item{a}{The value to test.}
}

\details{
Checks element-wise whether an object is equal to TRUE. If any of the elements are false, then an error will be thrown. If \sQuote{a} is TRUE, by default the function will print a statement informing the user the test has passed.  To turn off these print statements use \code{options("mxPrintUnitTests" = FALSE)}.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\seealso{\itemize{
\item{\code{\link{omxCheckCloseEnough}}}{for approximate equality testing within an epsilon.}
\item{\code{\link{omxCheckWithinPercentError}}}{for approximate equality testing within a percentage.}
\item{\code{\link{omxCheckEquals}}}{for exact equality testing.}
\item{\code{\link{omxCheckSetEquals}}}{for set equality testing.}
}}


\examples{

omxCheckTrue(1 + 1 == 2)

omxCheckTrue(matrix(TRUE, 3, 3))

# Throws an error
try(omxCheckTrue(FALSE))

}
