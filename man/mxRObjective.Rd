\name{mxRObjective}
\alias{mxRObjective}

\title{Function to Create MxRObjective Object}

\description{
   This function creates a new MxRObjective object.
}

\usage{
mxRObjective(objfun)
}

\arguments{
   \item{objfun}{A function that accepts two arguments.}
}

\details{
	The objfun argument must be a function that accepts two arguments. The first argument
	is the mxModel that should be evaluated, and the second argument is some persistent 
	state information that can be stored between one iteration of optimization to the next
 	iteration. It is valid for the function to simply ignore the second argument.

	The function must return either a single numeric value, or a list of exactly two elements.
	If the function returns a list, the first argument must be a single numeric value and the 
	second element will be the new persistent state information to be passed into this function
	at the next iteration. The single numeric value will be used by the optimizer to perform
	optimization.

	The initial default value for the persistant state information is NA.
}

\value{
Returns a new MxRObjective object.
}

\references{
The Mx manual can be found at http://www.vcu.edu/mx/documentation.html.
}

\examples{

A <- mxMatrix(nrow = 2, ncol = 2, values = c(1:4), free = TRUE, name = 'A')

squared <- function(x) { x ^ 2 }

objFunction <- function(model, state) {
	values <- model[['A']]@values 
	return(squared(values[1,1] - 4) + squared(values[1,2] - 3) +
		squared(values[2,1] - 2) + squared(values[2,2] - 1))
}
objective <- mxRObjective(objFunction)

model <- mxModel('model', A, objective)

}
