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
	name \tab - \tab The name of the object \cr
	min \tab - \tab The lower bound \cr 
	max \tab - \tab The upper bound \cr
	parameters \tab - \tab The vector of parameter names \cr     
  }
}
