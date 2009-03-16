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
}
