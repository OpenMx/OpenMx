\name{mxConstraint}
\alias{mxConstraint}

\title{Function To Create an MxConstraint Object}

\description{
This function creates a new \link{MxConstraint} object.
}

\usage{
mxConstraint(alg1, relation, alg2, name = NA)
}

\arguments{
\item{alg1}{A character string indicating the name of an \link{MxAlgebra} or \link{MxMatrix} object, whose relationship to the object specified in the 'alg2' argument is constrained.}
\item{relation}{A character string indicating the relation between 'alg1' and 'alg2'. Must be either "<", "=", or ">".}
\item{alg2}{A character string indicating the name of an \link{MxAlgebra} or \link{MxMatrix} object, whose relationship to the object specified in the 'alg1' argument is constrained.}   
\item{name}{An optional character string indicating the name of this object.}
}

\details{
The mxConstraint function defines relationships between two \link{MxAlgebra} or \link{MxMatrix} objects. They are used to affect the estimation of free parameters in the referenced objects. The relationships "<", "=" and ">" are supported. To affect an estimation or optimization, an \link{MxConstraint} object must be included in an \link{MxModel} object with all referenced \link{MxAlgebra} and \link{MxMatrix} objects.

The mxConstraint function may not be used to constrain free parameters, either by name or by their position in an \link{MxMatrix} or \link{MxAlgebra} object. Free parameters in the same \link{MxModel} may be constrained to equality by giving them the same name in their respective 'spec' matrices.
}

\value{
Returns an \link{MxConstraint} object. If used as an argument in an \link{MxModel} object, the objects referenced in the 'alg1' and 'alg2' arguments must also be included prior to optimization.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\seealso{
\link{MxConstraint} for the S4 class created by mxConstraint. \link{MxAlgebra} and \link{MxMatrix} for objects which may be entered as arguments in the 'alg' and 'alg2' arguments, and \link{mxAlgebra} and \link{mxMatrix} for the functions that create them. More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}

\examples{
#Create a constraint between MxMatrices 'A' and 'B'
constraint <- mxConstraint('A', '=', 'B', name = 'constraint')

#Constrain a matrix of free parameters 'F' to be greater than a matrix of zeros 'Z'
F <- mxMatrix(nrow=2, ncol=2, free=TRUE, name="F")
Z <- mxMatrix(type="Zero", nrow=2, ncol=2, name="Z")
positive <- mxConstraint("F", ">", "Z")
}
