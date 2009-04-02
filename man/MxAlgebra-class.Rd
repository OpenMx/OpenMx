\name{MxAlgebra-class}
\alias{MxAlgebra}

\title{MxAlgebra Class}

\description{
   MxAlgebra is an S4 class. An MxAlgebra object is a \link[=Named-entity]{named entity}.
   New instances of this class can be created using the function \link{mxAlgebra}.
}

\details{
The MxAlgebra class has the following slots:

  \tabular{rcl}{
	\tab \tab \cr
	name \tab - \tab The name of the object \cr
  formula \tab - \tab The R expression to be evaluated \cr 
	result \tab - \tab Either NULL or a 1x1 matrix with the computation result 
  }
  
The 'name' slot is the name of the MxAlgebra object. Use of MxAlgebra objects in the \link{mxConstraint} function or an \link{OpenMx} objective function requires reference by name.

The 'formula' slot is an expression containing one or more \link{MxMatrix} objects. These objects are operated on or related to one another using one or more operations detailed in the \link{mxAlgebra} help file.

The 'result' slot is used to hold the results of an optimization on the expression in the 'formula' slot. If this MxAlgebra has not been used as an argument in an objective function and subsequently included in an \link{MxModel} object and executed using the \link{mxRun} function, this slot will have a value of NULL. If this MxAlgebra has been used as an argument in an objective function and subsequently included in an \link{MxModel} object and executed using the \link{mxRun} function, then this slot will have a 1x1 matrix with the results of that optimization as its only value.

Slots may be referenced with the @ symbol. See the documentation for \link{Extract} and the examples in the \link{mxAlgebra} document for more information.
}

\references{
The Mx manual may be found at http://www.vcu.edu/mx/documentation.html.
}

\seealso{
\link{mxAlgebra} for creating MxAlgebra objects. \link{MxMatrix} and \link{mxMatrix} for objects which may be entered in the 'formula' slot and the function that creates them. More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}
