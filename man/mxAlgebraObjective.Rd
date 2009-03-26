\name{mxAlgebraObjective}
\alias{mxAlgebraObjective}

\title{Function to Create MxAlgebraObjective Object}

\description{
   This function creates a new MxAlgebraObjective object.
}

\usage{
mxAlgebraObjective(algebra, name = NA)
}

\arguments{
   \item{algebra}{A character string indicating the 
		name of the algebra or matrix to use for optimization.}
   \item{name}{A character string indicating the
   		name of this object. Optional.}
}

\value{
	Returns a new MxAlgebraObjective object.
}

\examples{

   A <- mxMatrix('Full', nrow = 3, ncol = 3, values = 1:9, name = 'A')

   B <- mxAlgebra(A + A, name = 'B')

   objective <- mxAlgebraObjective('B')

   model <- mxModel(A, B, objective)

}