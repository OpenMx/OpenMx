\name{mxFactor}
\alias{mxFactor}

\title{Fail-safe Factors}

\description{
   This is a wrapper for the R function \code{\link{factor}}. The \code{factor} function will skip missing levels {1,2,3,5} if they are not explicitly specified.  Use this function to force yourself to write out the levels.
}

\usage{
mxFactor(x = character(), levels, labels = levels, 
	exclude = NA, ordered = TRUE)
}

\arguments{
	\item{x}{a vector of data, usually taking a small number of distinct values.}
	\item{levels}{a mandatory vector of the values that 'x' might have taken.}
	\item{labels}{_either_ an optional vector of labels for the levels, _or_ a character string of length 1.}
	\item{exclude}{a vector of values to be excluded from the set of levels.}
	\item{ordered}{logical flag to determine if the levels should be regarded as ordered (in the order given).  Required to be TRUE.}
}

\details{
The function will throw an error if \sQuote{ordered} is not TRUE or if \sQuote{levels} is missing.  See \code{\link{factor}} for more information on creating ordered factors. 

}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{

ff <- mxFactor(substring("statistics", 1:10, 1:10), levels=letters)
as.integer(ff)  # the internal codes

factor(ff)      # drops the levels that do not occur.
                # AVOID DOING THIS UNINTENTIONALLY.
    
}
