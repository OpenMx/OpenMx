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
Objects can be created by calls of the form 
\code{mxExpectationGREML(V, yvars, Xvars, addOnes, blockByPheno, 
                        staggerZeroes, dataset.is.yX, casesToDropFromV)}.
}
\section{Slots}{
  \describe{
    \item{\code{V}:}{Object of class \code{"MxCharOrNumber"}.  Identifies the \code{\link{MxAlgebra}} or \code{\link{MxMatrix}} to serve as the 'V' matrix.}
    \item{\code{yvars}:}{Character vector.  Each string names a column of the raw dataset, to be used as a phenotypes.}
    \item{\code{Xvars}:}{A list of data column names, specifying the covariates to be used with each phenotype.}
    \item{\code{addOnes}:}{Logical; pertains to data-handling at runtime.}
    \item{\code{blockByPheno}:}{Logical; pertains to data-handling at runtime.}
    \item{\code{staggerZeroes}:}{Logical; pertains to data-handling at runtime.}
    \item{\code{dataset.is.yX}:}{Logical; pertains to data-handling at runtime.}
    \item{\code{y}:}{A matrix, to contain the 'y' vector of phenotypic observations.}
    \item{\code{X}:}{A matrix, to contain the 'X' matrix of covariates.}
    \item{\code{yXcolnames}:}{Character vector; used to store the column names of 'y' and 'X.'}
    \item{\code{casesToDrop}:}{Integer vector, specifying the rows and columns of the 'V' matrix to be removed at runtime.}
    \item{\code{b}:}{A matrix, to contain the vector of regression coefficients calculated at runtime.}
    \item{\code{bcov}:}{A matrix, to contain the sampling covariance matrix of the regression coefficients calculated at runtime.}
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
See \code{\link{mxExpectationGREML}()} for creating MxExpectationGREML objects, and for more information generally concerning GREML analyses, including a complete example. More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}
\examples{
showClass("MxExpectationGREML")
}
\keyword{classes}
