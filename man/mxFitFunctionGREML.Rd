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
   augGrad=character(0), augHess=character(0), 
   autoDerivType=c("semiAnalyt","numeric"),infoMatType=c("average","expected"))
}

\arguments{
  \item{dV}{Vector of character strings; defaults to a character vector of length zero.  If a value of non-zero length is provided, it must be a \emph{named} character vector.  This vector's names must be the labels of each free parameter in the model.  The vector's elements (i.e., the character strings themselves) must be the names of \code{\link{MxAlgebra}} or \code{\link{MxMatrix}} object(s), each of which equals the first partial derivative of the 'V' matrix with respect to the corresponding free parameter.}
  \item{aug}{Character string; defaults to a character vector of length zero.  Any elements after the first are ignored.  The string should name a 1x1 \code{\link{MxMatrix}} or an \code{\link{MxAlgebra}} that evaluates to a 1x1 matrix.  The named object will be used as an "augmentation" to the GREML fitfunction--specifically, the [1,1] value of the object named by \code{aug} will be added to the GREML fitfunction value at each function evaluation during optimization.  The augmentation can be used to regularize estimation with a prior likelihood, or to use penalty functions to approximate constraints.}
  \item{augGrad}{Character string; defaults to a character vector of length zero.  Any elements after the first are ignored.  The string should name a \code{\link{MxMatrix}} or an \code{\link{MxAlgebra}} that evaluates to the gradient of \code{aug} with respect to free parameters.  The gradient can be either a column or row vector.  The free parameters corresponding to the elements of the gradient vector are taken from the names of argument \code{dV}, e.g. if the third name of \code{dV} is \code{'va'}, then the third element of the gradient vector should be the first partial derivative of the augmentation function with respect to \code{'va'}.  Ignored unless both \code{dV} and \code{aug} have nonzero length.}
  \item{augHess}{Character string; defaults to a character vector of length zero.  Any elements after the first are ignored.  The string should name a \code{\link{MxMatrix}} or an \code{\link{MxAlgebra}} that evaluates to the Hessian of \code{aug} with respect to free parameters.  The free parameters corresponding to each row and column of this matrix are dictated by the names of argument \code{dV}, in the same manner as for the elements of \code{augGrad}.  Ignored unless both \code{dV} and \code{aug} have nonzero length.  Providing a nonzero-length value for \code{augHess} but not \code{augGrad} will result in an error at runtime.}
  \item{autoDerivType}{"Automatic derivative type."  Character string, either "semiAnalyt" (default) or "numeric".  See details below.}
  \item{infoMatType}{"Information matrix type."  Character string, either "average" (default) or "expected".  See details below.}
}
\details{
Making effective use of arguments \code{dV}, \code{augGrad}, and \code{augHess} may require a custom \code{\link{mxComputeSequence}()}.  The derivatives of the REML loglikelihood function with respect to parameters can be internally computed from the derivatives of the 'V' matrix supplied via \code{dV}.  The loglikelihood's first derivatives thus computed will always be exact, but its matrix of second partial derivatives (i.e., its Hessian matrix) will be approximated by either the average or expected information matrix, per the value of argument \code{infoMatType}.  The average information matrix is faster to compute, but may not provide a good approximation to the Hessian if 'V' is not linear in the model's free parameters.  The expected information matrix is slower to compute, but does not assume that 'V' is linear in the free parameters.  Neither information matrix will be a good approximation to the Hessian unless the derivatives of 'V' evaluate to symmetric matrices the same size as 'V'.  Note also that these loglikelihood derivatives do not reflect the influence of any parameter bounds or \code{\link{MxConstraint}}s. Internally, the derivatives of the 'V' matrix are assumed to be symmetric, and the elements above their main diagonals are ignored.

Formerly, if any derivatives were provided via \code{dV}, then derivatives had to be provided for \emph{every} free parameter in the MxModel.  Currently, users may provide derivatives of 'V' via \code{dV} with respect to some or all free parameters.  Note that the gradient and Hessian of the augmentation must be complete, i.e. contain derivatives of the augmentation with respect to every parameter or pair of parameters respectively.

If there are any free parameters with respect to which the user did not provide an analytic derivative of 'V', OpenMx will automatically calculate the necessary loglikelihood derivatives according to \code{autoDerivType}.  If \code{autoDerivType="semiAnalyt"}, the GREML fitfunction backend will calculate the missing derivatives in a "semi-analytic" fashion.  Specifically, the backend will numerically differentiate 'V' with respect to the relevant parameter(s), and use those numeric matrix derivatives to analytically calculate the needed loglikelihood derivatives.  If \code{autoDerivType="numeric"}, the needed loglikelihood derivatives will be calculated numerically, via finite-differences.

Argument \code{aug} is intended to allow users to provide penalty functions or prior likelihoods in order to approximate constraints or to regularize optimization.  The user is warned that careless use of this augmentation feature may undermine the validity of his/her statistical inferences.
}
\value{
Returns a new object of class \code{\link{MxFitFunctionGREML}}.
}

\references{
Kirkpatrick RM, Pritikin JN, Hunter MD, & Neale MC.  (2021).  Combining structural-equation modeling with genomic-relatedness matrix restricted maximum likelihood in OpenMx.  Behavior Genetics 51: 331-342.  \url{https://doi.org/10.1007/s10519-020-10037-5}

The first software implementation of "GREML":\cr
Yang J, Lee SH, Goddard ME, Visscher PM.  (2011).  GCTA: a tool for genome-wide complex trait analysis.  American Journal of Human Genetics 88: 76-82.  \url{https://doi.org/10.1016/j.ajhg.2010.11.011}

One of the first uses of the acronym "GREML":\cr
Benjamin DJ, Cesarini D, van der Loos MJHM, Dawes CT, Koellinger PD, et al. (2012).  The genetic architecture of economic and political preferences. Proceedings of the National Academy of Sciences 109: 8026-8031. doi: 10.1073/pnas.1120666109

The OpenMx User's guide can be found at \url{https://openmx.ssri.psu.edu/documentation/}.
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
