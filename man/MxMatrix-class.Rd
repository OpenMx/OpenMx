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
	}