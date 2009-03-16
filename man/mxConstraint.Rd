\name{mxConstraint}
\alias{mxConstraint}

\title{Function To Create MxConstraint Object}

\description{
   This function creates a new \link{MxConstraint} object.
}

\usage{
	mxConstraint(alg1, relation, alg2, name = NA)
}


\arguments{
   \item{alg1}{A character string indicating the 
   		name of an \link{MxAlgebra} or \link{MxMatrix} object.}
   \item{relation}{A character string indicating the relation between
   		alg1 and alg2. Must be either "<", "=", or ">".}
   \item{alg2}{A character string indicating the
   		name of an \link{MxAlgebra} or \link{MxMatrix} object.}   
   \item{name}{A character string indicating the
   		name of this object. Optional.}
}

\value{
	Returns a new \link{MxConstraint} object.
}

\examples{

   #Create a constraint between MxMatrices 'A' and 'B'
   constraint <- mxConstraint('A', '=', 'B', name = 'constraint')

}
