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
	result \tab - \tab a matrix with the computation result
  }
  
The \sQuote{name} slot is the name of the MxAlgebra object. Use of MxAlgebra objects in the \code{\link{mxConstraint}} function or an objective function requires reference by name.

The \sQuote{formula} slot is an expression containing one or more \code{\link{MxMatrix}} objects. These objects are operated on or related to one another using one or more operations detailed in the \code{\link{mxAlgebra}} help file.

The \sQuote{result} slot is used to hold the results of computing the expression in the \sQuote{formula} slot. If the containing model has not been executed, then the \sQuote{result} slot will hold a 0 x 0 matrix.  Otherwise the slot will store the computed value of the algebra using the final estimates of the free parameters.

Slots may be referenced with the @ symbol. See the documentation for \link[methods]{Classes} and the examples in the \code{\link{mxAlgebra}} document for more information.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.}

\seealso{
\code{\link{mxAlgebra}}, \code{\link{mxMatrix}}, \code{\link{MxMatrix}} 
}
