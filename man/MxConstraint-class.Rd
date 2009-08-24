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
	alg1 \tab - \tab The name of an \link{MxAlgebra} or \link{MxMatrix} object \cr 
	relation \tab - \tab A character string, either '<', '=', or '>' \cr
	alg2 \tab - \tab The name of an \link{MxAlgebra} or \link{MxMatrix} object \cr     
  }
  
The 'name' slot is the name of the MxConstraint object. %Use of MxConstraint objects in other functions in the \link{OpenMx} library may require reference by name.

The 'alg1' and 'alg2' slots hold \link{MxAlgebra} or \link{MxMatrix} objects whose relationship is constrained or defined by the contents of the 'relation' slot. To affect an estimation or optimization, an MxConstraint object must be included in an \link{MxModel} object with all referenced \link{MxAlgebra} and \link{MxMatrix} objects.

Slots may be referenced with the @ symbol. See the documentation for \link[methods]{Classes} and the examples in the \link{mxConstraint} document for more information.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\seealso{
\link{mxConstraint} for the function that creates MxConstraint objects. \link{MxAlgebra} and \link{MxMatrix} for objects which may be entered as arguments in the 'alg' and 'alg2' arguments, and \link{mxAlgebra} and \link{mxMatrix} for the functions that create them. More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}
