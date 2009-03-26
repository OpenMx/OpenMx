\name{mxBounds}
\alias{mxBounds}

\title{Function To Create MxBounds Object}

\description{
   This function creates a new \link{MxBounds} object.
}

\usage{
mxBounds(parameters, min = NA, max = NA, name = NA)
}

\arguments{
   \item{parameters}{A character vector indicating the 
   		names of the parameters on which to apply bounds.}
   \item{min}{A numeric value for the lower bound. NA means use default value.}
   \item{max}{A numeric value for the upper bound. NA means use default value.}   
   \item{name}{A character string indicating the
   		name of this object. Optional.}
}

\value{
	Returns a new \link{MxBounds} object.
}

\examples{

   #Create a bound for parameters 'A' and 'B'
   bounds <- mxBounds(c('A', 'B'), 3, 5, name = 'ABbounds')

}
