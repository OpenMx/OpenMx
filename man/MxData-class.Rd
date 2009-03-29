\name{MxData-class}
\alias{MxData}

\title{MxData Class}

\description{
   MxData is an S4 class. An MxData object is a \link[=Named-entity]{named entity}.
   New instances of this class can be created using the function \link{mxData}.
   MxData is an S4 class union. An MxData object is either \link{NULL} or a 
   MxNonNullData object.
}

\details{
  The MxNonNullData class has the following slots:

  \tabular{rcl}{
	\tab \tab \cr
	name \tab - \tab The name of the object \cr
	matrix \tab - \tab Either a matrix or a data frame \cr
	vector \tab - \tab A vector for means, or NA if missing \cr
	type \tab - \tab Either 'raw', 'cov', 'cor', or 'sscp' \cr
	numObs \tab - \tab The number of oberservations \cr
  }
}
