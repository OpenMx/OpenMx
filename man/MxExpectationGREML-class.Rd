\name{MxExpectationGREML-class}
\Rdversion{1.1}
\docType{class}
\alias{MxExpectationGREML-class}
\alias{MxExpectationGREML}

\title{Class "MxExpectationGREML"}
\description{
\code{MxExpectationGREML} is a type of expectation class.  It contains the necessary elements for specifying a GREML model.  For more information, see \code{\link{mxExpectationGREML}()}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{mxExpectationGREML(V, X, y, dV, fixedEffects)}.
}
\section{Slots}{
  \describe{
    \item{\code{V}:}{Object of class \code{"MxCharOrNumber"}.  Identifies the \code{\link{MxAlgebra}} or \code{\link{MxMatrix}} to serve as the 'V' matrix.}
    \item{\code{X}:}{Object of class \code{"MxCharOrNumber"}.  Identifies the \code{\link{MxAlgebra}} or \code{\link{MxMatrix}} to serve as the 'X' matrix.}
    \item{\code{y}:}{Object of class \code{"MxCharOrNumber"}.  Identifies the \code{\link{MxAlgebra}} or \code{\link{MxMatrix}} to serve as the 'y' vector.}
    \item{\code{dV}:}{Object of class \code{"MxCharOrNumber"}.  Identifies the \code{\link{MxAlgebra}} or \code{\link{MxMatrix}} object(s) to serve as the derivatives of 'V' with respect to free parameters.}
    \item{\code{dVnames}:}{Vector of character strings; names of the free parameters corresponding to slot \code{dV}}
    \item{\code{numFixEff}:}{Integer number of covariates in 'X.'}
    \item{\code{dims}:}{Object of class \code{"character"}.}
    \item{\code{definitionVars}:}{Object of class \code{"list"}.}
    \item{\code{numStats}:}{Numeric; number of observed statistics.}
    \item{\code{name}:}{Object of class \code{"character"}.}
    \item{\code{data}:}{Object of class \code{"MxCharOrNumber"}.}
    \item{\code{submodels}:}{Object of class \code{"MxOptionalCharOrNumber"}.}
    \item{\code{container}:}{Object of class \code{"MxOptionalCharOrNumber"}.}
    \item{\code{.runDims}:}{Object of class \code{"character"}.}
  }
}
\section{Extends}{
Class \code{"MxBaseExpectation"}, directly.
Class \code{"MxBaseNamed"}, by class "MxBaseExpectation", distance 2.
Class \code{"MxExpectation"}, by class "MxBaseExpectation", distance 2.
}
\section{Methods}{
No methods defined with class "MxExpectationGREML" in the signature.
}
\references{
The OpenMx User's guide can be found at \url{http://openmx.psyc.virginia.edu/documentation}.
}

\seealso{
See \code{\link{mxExpectationGREML}()} for creating MxExpectationGREML objects, and for more information generally concerning GREML analyses. More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}
\examples{
showClass("MxExpectationGREML")
}
\keyword{classes}
