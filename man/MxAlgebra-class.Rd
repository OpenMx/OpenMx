\name{MxAlgebra-class}
\alias{MxAlgebra-class}
\alias{MxAlgebra}

\title{MxAlgebra Class}

\description{
   MxAlgebra is an S4 class. An MxAlgebra object is a named entity.
   New instances of this class can be created using the function \link{mxAlgebra}.
}

\details{
  The MxAlgebra class has the following slots:
  \tabular{rcl}{
	name \tab - \tab The name of the object \cr
    formula \tab - \tab The R expression to be evaluated \cr 
	result \tab - \tab Either NULL or a 1x1 matrix with the computation result 
  }
}
