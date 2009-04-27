\name{MxBounds-class}
\alias{MxBounds}

\title{MxBounds Class}

\description{
   MxBounds is an S4 class. An MxBounds object 
   is a \link[=Named-entity]{named entity}. New instances of this class can
   be created using the function \link{mxBounds}.
}

\details{
  The MxBounds class has the following slots:

  \tabular{rcl}{
	\tab \tab \cr
	min \tab - \tab The lower bound \cr 
	max \tab - \tab The upper bound \cr
	parameters \tab - \tab The vector of parameter names \cr     
  }
  
The 'min' and 'max' slots hold scalar numeric values for the lower and upper bounds on the list of parameters, respectively. If a vector of values is input, only the first value is used.

Parameters may be any free parameter or parameters from an \link{MxMatrix} object. Parameters may be referenced either by name or by referring to their position in the 'spec' matrix of an \code{MxMatrix} object. To affect an estimation or optimization, an MxBounds object must be included in an \link{MxModel} object with all referenced \link{MxAlgebra} and \link{MxMatrix} objects.

Slots may be referenced with the @ symbol. See the documentation for \link{Extract} and the examples in the \link{mxBounds} document for more information.
}

\references{
More information about Mx may be found at http://www.vcu.edu/mx/documentation.html.
}

\seealso{
\link{mxBounds} for the function that creates MxBounds objects. \link{MxMatrix} and \link{mxMatrix} for free parameter specification. More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}
