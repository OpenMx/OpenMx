\name{mxFitFunctionGREML}
\alias{mxFitFunctionGREML}

\title{
Create MxFitFunctionGREML Object
}
\description{
This function creates a new \code{\link{MxFitFunctionGREML}} object.
}
\usage{
mxFitFunctionGREML(dV=character(0))
}

\arguments{
  \item{dV}{Vector of character strings; defaults to a character vector of length zero.  If a value of non-zero length is provided, it must be a \emph{named} character vector.  This vector's names must be the labels of free parameters in the model.  The vector's elements (i.e., the character strings themselves) must be the names of \code{\link{MxAlgebra}} or \code{\link{MxMatrix}} object(s), each of which equals the first partial derivative of the 'V' matrix with respect to the corresponding free parameter.}
}
\details{
Making effective use of argument \code{dV} will usually require a custom \code{\link{mxComputeSequence}()}.  The derivatives of the REML loglikelihood function with respect to parameters can be internally computed from the derivatives of the 'V' matrix supplied via \code{dV}.  These loglikelihood derivatives will be valid as long as (1) the derivatives of 'V' evaluate to symmetric matrices the same size as 'V,' and (2) the model contains no \code{\link{MxConstraint}}s. Internally, the derivatives of the 'V' matrix are assumed to be symmetric, and the elements above their main diagonals are ignored.
}
\value{
Returns a new object of class \code{\link{MxFitFunctionGREML}}.
}
\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\seealso{
See \code{\link{MxFitFunctionGREML}} for the S4 class created by \code{mxFitFunctionGREML()}.  For more information generally concerning GREML analyses, including a complete example, see \code{\link{mxExpectationGREML}()}.  More information about the OpenMx package may be found \link[=OpenMx]{here}.
}
\examples{
gff <- mxFitFunctionGREML()
str(gff)
}
