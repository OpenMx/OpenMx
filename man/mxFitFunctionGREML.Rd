\name{mxFitFunctionGREML}
\alias{mxFitFunctionGREML}

\title{
Create MxFitFunctionGREML Object
}
\description{
This function creates a new \code{\link{MxFitFunctionGREML}} object.
}
\usage{
mxFitFunctionGREML(dV=character(0), aug=character(0),
   augGrad=character(0), augHess=character(0))
}

\arguments{
  \item{dV}{Vector of character strings; defaults to a character vector of length zero.  If a value of non-zero length is provided, it must be a \emph{named} character vector.  This vector's names must be the labels of each free parameter in the model.  The vector's elements (i.e., the character strings themselves) must be the names of \code{\link{MxAlgebra}} or \code{\link{MxMatrix}} object(s), each of which equals the first partial derivative of the 'V' matrix with respect to the corresponding free parameter.}
  \item{aug}{Character string; defaults to a character vector of length zero.  Any elements after the first are ignored.  The string should name a 1x1 \code{\link{MxMatrix}} or an \code{\link{MxAlgebra}} that evaluates to a 1x1 matrix.  The named object will be used as an "augmentation" to the GREML fitfunction--specifically, the [1,1] value of the object named by \code{aug} will be added to the GREML fitfunction value at each function evaluation during optimization.  The augmentation can be used to regularize estimation with a prior likelihood, or to use penalty functions to approximate constraints.}
  \item{augGrad}{Character string; defaults to a character vector of length zero.  Any elements after the first are ignored.  The string should name a \code{\link{MxMatrix}} or an \code{\link{MxAlgebra}} that evaluates to the gradient of \code{aug} with respect to free parameters.  The gradient can be either a column or row vector.  The free parameters corresponding to the elements of the gradient vector are taken from the names of argument \code{dV}, e.g. if the third name of \code{dV} is \code{'va'}, then the third element of the gradient vector should be the first partial derivative of the augmentation function with respect to \code{'va'}.  Ignored unless both \code{dV} and \code{aug} have nonzero length.}
  \item{augHess}{Character string; defaults to a character vector of length zero.  Any elements after the first are ignored.  The string should name a \code{\link{MxMatrix}} or an \code{\link{MxAlgebra}} that evaluates to the Hessian of \code{aug} with respect to free parameters.  The free parameters corresponding to each row and column of this matrix are dictated by the names of argument \code{dV}, in the same manner as for the elements of \code{augGrad}.  Ignored unless both \code{dV} and \code{aug} have nonzero length.  Providing a nonzero-length value for \code{augHess} but not \code{augGrad} will result in an error at runtime.}
}
\details{
Making effective use of arguments \code{dV}, \code{augGrad}, and \code{augHess} will usually require a custom \code{\link{mxComputeSequence}()}.  The derivatives of the REML loglikelihood function with respect to parameters can be internally computed from the derivatives of the 'V' matrix supplied via \code{dV}.  These loglikelihood derivatives will be valid as long as (1) the derivatives of 'V' evaluate to symmetric matrices the same size as 'V,' and (2) the model contains no \code{\link{MxConstraint}}s. Internally, the derivatives of the 'V' matrix are assumed to be symmetric, and the elements above their main diagonals are ignored.

Currently, if any derivatives are provided via \code{dV}, then derivatives must be provided for \emph{every} free parameter in the MxModel.

Argument \code{aug} is intended to allow users to provide penalty functions or prior likelihoods in order to approximate constraints or to regularize optimization.  The user is warned that careless use of this augmentation feature may undermine the validity of his/her statistical inferences.
}
\value{
Returns a new object of class \code{\link{MxFitFunctionGREML}}.
}
\references{
The OpenMx User's guide can be found at http://openmx.ssri.psu.edu/documentation.
}

\seealso{
See \code{\link{MxFitFunctionGREML}} for the S4 class created by \code{mxFitFunctionGREML()}.  For more information generally concerning GREML analyses, including a complete example, see \code{\link{mxExpectationGREML}()}.

Other fit functions:
\code{\link{mxFitFunctionMultigroup}}, \code{\link{mxFitFunctionML}}, 
\code{\link{mxFitFunctionWLS}}, \code{\link{mxFitFunctionAlgebra}}, 
\code{\link{mxFitFunctionR}}, \code{\link{mxFitFunctionRow}}

More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}

\examples{
gff <- mxFitFunctionGREML()
str(gff)
}
