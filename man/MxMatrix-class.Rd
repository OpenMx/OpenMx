\name{MxMatrix-class}
\alias{MxMatrix}
	
\title{MxMatrix Class}

\description{
   MxMatrix is an S4 class. An MxMatrix object is a \link[=Named-entity]{named entity}.
   New instances of this class can be created using the function \link{mxMatrix}.
   MxMatrix objects may be used as arguments in other functions from the OpenMx library, including \link{mxAlgebra}, \link{mxConstraint}, and \link{mxModel}.
}
	
	\details{
The MxMatrix class has the following slots:
	  \tabular{rcl}{
	        \tab \tab \cr
	        name \tab - \tab The name of the object \cr
	        specification \tab - \tab  The specification matrix, \cr
	        values \tab - \tab The values matrix
	  }
	  
The 'name' slot is the name of the MxMatrix object. Use of MxMatrix objects in an \link{mxAlgebra} or \link{mxConstraint} function requires reference by name.

The 'specification' slot takes a matrix which describes the location of free and fixed parameters. Fixed parameters are specified with numeric zeros, while free parameters may be specified with character strings, non-zero numeric values, and 'NA's. Free parameters are elements of an MxMatrix object whose values may be changed by an objective function when that MxMatrix object is included in an \link{MxModel} object and evaluated using the \link{mxRun} function.

The 'values' slot takes a matrix of numeric values. If an element is specified as a fixed parameter in the 'specification' matrix, then the element in the 'values' matrix is treated as a constant value and cannot be altered or updated by an objective function when included in an \link{mxRun} function. If an element is specified as a free parameter in the 'specification' matrix, the element in the 'value' matrix is considered a starting value and can be changed by an objective function when included in an \link{mxRun} function. 
}

\references{
The Mx manual may be found at http://www.vcu.edu/mx/documentation.html.
}

\seealso{
\link{mxMatrix} for creating MxMatrix objects. More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}