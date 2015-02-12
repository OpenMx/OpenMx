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
Objects can be created by calls of the form \code{mxFitFunctionGREML(casesToDrop, dropNAfromV)}.
}
\section{Slots}{
  \describe{
    \item{\code{casesToDrop}:}{Integer; the indices of the rows of 'y' and 'X' that were dropped because of missing values.}
    \item{\code{dropNAfromV}:}{Logical; should the rows (and columns) indexed in slot \code{casesToDrop} also be internally dropped from 'V' at runtime?}
    \item{\code{info}:}{Object of class \code{"list"}.  As of this writing, this slot contains the output for the fixed effects.}
    \item{\code{dependencies}:}{Object of class \code{"integer"}.}
    \item{\code{expectation}:}{Object of class \code{"integer"}.}
    \item{\code{vector}:}{Object of class \code{"logical"}.}
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
The OpenMx User's guide can be found at \url{http://openmx.psyc.virginia.edu/documentation}.
}

\seealso{
See \code{\link{mxFitFunctionGREML}()} for creating MxFitFunctionGREML objects.  See \code{\link{mxExpectationGREML}()} for creating MxExpectationGREML objects, and for more information generally concerning GREML analyses. More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}
\examples{
showClass("MxFitFunctionGREML")
}
\keyword{classes}
