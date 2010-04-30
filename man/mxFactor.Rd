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
	\item{x}{either a vector of data or a data.frame object.}
	\item{levels}{a mandatory vector of the values that 'x' might have taken.}
	\item{labels}{_either_ an optional vector of labels for the levels, _or_ a character string of length 1.}
	\item{exclude}{a vector of values to be excluded from the set of levels.}
	\item{ordered}{logical flag to determine if the levels should be regarded as ordered (in the order given).  Required to be TRUE.}
}

\details{
If \sQuote{x} is a data.frame, then all of the columns of \sQuote{x} are converted into ordered factors. If \sQuote{x} is a data.frame, then \sQuote{levels} and \sQuote{labels} may be either a list or a vector. When \sQuote{levels} is a list, then different levels are assigned to different columns of the constructed data.frame object. When \sQuote{levels} is a vector, then the same levels are assigned to all the columns of the data.frame object. The function will throw an error if \sQuote{ordered} is not TRUE or if \sQuote{levels} is missing.  See \code{\link{factor}} for more information on creating ordered factors. 

}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{

ff <- mxFactor(substring("statistics", 1:10, 1:10), levels=letters)
as.integer(ff)  # the internal codes

factor(ff)      # drops the levels that do not occur.
                # AVOID DOING THIS UNINTENTIONALLY.

foo <- data.frame(x=c(1:3),y=c(4:6),z=c(7:9))
mxFactor(foo, c(1:9)) # same levels applied to all columns
mxFactor(foo, list(c(1:3), c(4:6), c(7:9))) # different levels to different columns

    
}
