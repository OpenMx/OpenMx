\name{MxConstraint-class}
\alias{MxConstraint}

\title{MxConstraint Class}

\description{
   MxConstraint is an S4 class. An MxConstraint object 
   is a \link[=Named-entity]{named entity}. New instances of this class can
   be created using the function \link{mxConstraint}.
}

\details{
The MxConstraint class has the following slots:

  \tabular{rcl}{
	\tab \tab \cr
	name \tab - \tab The name of the object \cr
	formula \tab - \tab The R expression to be evaluated
  }
  
The \sQuote{name} slot is the name of the MxConstraint object. Use of MxConstraint objects in other functions in the \link{OpenMx} library may require reference by name.

The \sQuote{formula} slot is an expression containing the expression to be evaluated. These objects are operated on or related to one another using one or more operations detailed in the \code{\link{mxConstraint}} help file.

Slots may be referenced with the @ symbol. See the documentation for \link[methods]{Classes} and the examples in the \link{mxConstraint} document for more information.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\seealso{
\link{mxConstraint} for the function that creates MxConstraint objects. \link{MxAlgebra} and \link{MxMatrix} for objects which may be entered as arguments in the 'alg' and 'alg2' arguments, and \link{mxAlgebra} and \link{mxMatrix} for the functions that create them. More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}
