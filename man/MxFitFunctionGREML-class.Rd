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
    \item{\code{dVnames}:}{Vector of character strings; names of the free parameters corresponding to slot \code{dV}.}
    \item{\code{info}:}{Object of class \code{"list"}.}
    \item{\code{dependencies}:}{Object of class \code{"integer"}.}
    \item{\code{expectation}:}{Object of class \code{"integer"}.}
    \item{\code{vector}:}{Object of class \code{"logical"}.}
    \item{\code{result}:}{Object of class \code{"matrix"}.}
    \item{\code{name}:}{Object of class \code{"character"}.}
    \item{\code{numObs}:}{Object of class \code{"integer"}.}
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
The OpenMx User's guide can be found at \url{http://openmx.psyc.virginia.edu/documentation}.
}

\seealso{
See \code{\link{mxFitFunctionGREML}()} for creating MxFitFunctionGREML objects.  See \code{\link{mxExpectationGREML}()} for creating MxExpectationGREML objects, and for more information generally concerning GREML analyses, including a complete example. More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}
\examples{
showClass("MxFitFunctionGREML")
}
\keyword{classes}
