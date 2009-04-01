\name{mxBounds}
\alias{mxBounds}

\title{Function To Create an \link{MxBounds} Object}

\description{
This function creates a new \link{MxBounds} object.
}

\usage{
mxBounds(parameters, min = NA, max = NA, name = NA)
}

\arguments{
\item{parameters}{A character vector indicating the names of the parameters on which to apply bounds.}
\item{min}{A numeric value for the lower bound. NA means use default value.}
\item{max}{A numeric value for the upper bound. NA means use default value.}   
\item{name}{An optional character string indicating the name of this object.}
}

\details{
Creates a set of boundaries or limits for a parameter or set of parameters. Parameters may be any free parameter or parameters from an \link{MxMatrix} object. Parameters may be referenced either by name or by referring to their position in the 'specification' matrix of an \code{MxMatrix} object.

Minima and maxima may be specified as scalar numeric values. If a vector of values is input, only the first value is used.
}
	
\value{
Returns a new \link{MxBounds} object. If used as an argument in an \link{MxModel} object, the parameters referenced in the 'parameters' argument must also be included prior to optimization.
}

\references{
More information about Mx may be found at http://www.vcu.edu/mx/documentation.html
}

\examples{
#Create lower and upper bounds for parameters 'A' and 'B'
bounds <- mxBounds(c('A', 'B'), 3, 5, name = 'ABbounds')

#Create a lower bound of zero for a set of variance parameters
varianceBounds <- mxBounds(c('Var1', 'Var2', 'Var3'), 0,
                           name = 'varianceBounds')
}