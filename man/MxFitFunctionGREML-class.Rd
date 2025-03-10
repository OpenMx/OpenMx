\name{MxFitFunctionGREML-class}
\Rdversion{1.1}
\docType{class}
\alias{MxFitFunctionGREML-class}
\alias{MxFitFunctionGREML}

\title{Class \code{"MxFitFunctionGREML"}}
\description{
\code{MxFitFunctionGREML} is the fitfunction class for GREML analyses.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{mxFitFunctionGREML(dV)}.
}
\section{Slots}{
  \describe{
    \item{\code{dV}:}{Object of class \code{"MxCharOrNumber"}.  Identifies the \code{\link{MxAlgebra}} or \code{\link{MxMatrix}} object(s) to serve as the derivatives of 'V' with respect to free parameters.}
    \item{\code{dVnames}:}{Vector of character strings; names of the free parameters corresponding to slot \code{dV}.  Provides the mapping from elements of \code{augGrad} and of rows and columns of \code{augHess} to free parameters.}
    \item{\code{MLfit}:}{Object of class \code{"numeric"}, equal to the maximum-likelihood fitfunction value (as opposed to the restricted maximum-likelihood value).}
    \item{\code{numObsAdjust}:}{Object of class \code{"integer"}.
      Number of observations adjustment.}
    \item{\code{aug}:}{Object of class \code{"MxCharOrNumber"}.  Identifies the \code{\link{MxAlgebra}} or \code{\link{MxMatrix}} object used to "augment" the fitfunction value at each function evaluation during optimization.}
    \item{\code{augGrad}:}{Object of class \code{"MxCharOrNumber"}.  Identifies the \code{\link{MxAlgebra}} or \code{\link{MxMatrix}} object(s) to serve as the first derivatives of \code{aug} with respect to free parameters.}
    \item{\code{augHess}:}{Object of class \code{"MxCharOrNumber"}.  Identifies the \code{\link{MxAlgebra}} or \code{\link{MxMatrix}} object(s) to serve as the second derivatives of \code{aug} with respect to free parameters.}
    \item{\code{autoDerivType}:}{Object of class \code{"character"}.  Dictates whether fitfunction derivatives automatically calculated by OpenMx should be numerical or "semi-analytic."}
    \item{\code{infoMatType}:}{Object of class \code{"character"}.  Dictates whether to calculate the average- or expected-information matrix.}
    \item{\code{dyhat}:}{Object of class \code{"MxCharOrNumber"}.  Identifies the \code{\link{MxAlgebra}} or \code{\link{MxMatrix}} object(s) to serve as the derivatives of 'yhat' with respect to free parameters.}
    \item{\code{dyhatnames}:}{Vector of character strings; names of the free parameters corresponding to slot \code{dyhat}.}
    \item{\code{dNames}:}{Vector of character strings.  If \code{dV} and/or \code{dyhat} have nonzero length, then populated at runtime by the set of non-redundant free-parameter names appearing in \code{dVnames} and/or \code{dyhatnames}.}
    \item{\code{.parallelDerivScheme}:}{Object of class \code{"integer"}.  Dictates how computation of the information matrix will be distributed across multiple threads.  This slot is intended for use in testing and debugging by the OpenMx developers, and not to be modified by the end user.}
    \item{\code{info}:}{Object of class \code{"list"}.}
    \item{\code{dependencies}:}{Object of class \code{"integer"}.}
    \item{\code{expectation}:}{Object of class \code{"integer"}.}
    \item{\code{vector}:}{Object of class \code{"logical"}.}
    \item{\code{rowDiagnostics}:}{Object of class \code{"logical"}.}
    \item{\code{result}:}{Object of class \code{"matrix"}.}
    \item{\code{name}:}{Object of class \code{"character"}.}

  }
}
\section{Extends}{
Class \code{"MxBaseFitFunction"}, directly.
Class \code{"MxBaseNamed"}, by class "MxBaseFitFunction", distance 2.
Class \code{"MxFitFunction"}, by class "MxBaseFitFunction", distance 2.
}
\section{Methods}{
No methods defined with class "MxFitFunctionGREML" in the signature.
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
See \code{\link{mxFitFunctionGREML}()} for creating MxFitFunctionGREML objects.  See \code{\link{mxExpectationGREML}()} for creating MxExpectationGREML objects, and for more information generally concerning GREML analyses, including a complete example. More information about the OpenMx package may be found \link[=OpenMx]{here}.
}
\examples{
showClass("MxFitFunctionGREML")
}
\keyword{classes}
