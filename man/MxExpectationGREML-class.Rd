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
                        staggerZeroes, dataset.is.yX, casesToDropFromV,
                        REML, yhat)}.
}
\section{Slots}{
  \describe{
    \item{\code{V}:}{Object of class \code{"MxCharOrNumber"}.  Identifies the \code{\link{MxAlgebra}} or \code{\link{MxMatrix}} to serve as the 'V' matrix.}
    \item{\code{yvars}:}{Character vector.  Each string names a column of the raw dataset, to be used as a phenotype.}
    \item{\code{Xvars}:}{A list of data column names, specifying the covariates to be used with each phenotype.}
    \item{\code{addOnes}:}{Logical; pertains to data-handling at runtime.}
    \item{\code{blockByPheno}:}{Logical; pertains to data-handling at runtime.}
    \item{\code{staggerZeroes}:}{Logical; pertains to data-handling at runtime.}
    \item{\code{dataset.is.yX}:}{Logical; pertains to data-handling at runtime.}
    \item{\code{y}:}{Object of class \code{"MxData"}.  Its \code{observed} slot will contain the phenotype vector, 'y.'}
    \item{\code{X}:}{A matrix, to contain the 'X' matrix of covariates.}
    \item{\code{yXcolnames}:}{Character vector; used to store the column names of 'y' and 'X.'}
    \item{\code{casesToDrop}:}{Integer vector, specifying the rows and columns of the 'V' matrix to be removed at runtime.}
    \item{\code{b}:}{A matrix, to contain the vector of regression coefficients calculated at runtime.}
    \item{\code{bcov}:}{A matrix, to contain the sampling covariance matrix of the regression coefficients calculated at runtime.}
    \item{\code{numFixEff}:}{Integer number of covariates in 'X.'}
    \item{\code{REML}:}{Logical.  Should restricted maximum-likelihood estimation be used?}
    \item{\code{yhat}:}{Object of class \code{"MxCharOrNumber"}.  Identifies the \code{\link{MxAlgebra}} or \code{\link{MxMatrix}} to serve as the model-expected phenotypic mean vector.}
    \item{\code{fitted.values}:}{A matrix, specifically, the column vector of model-predicted phenotypic means calculated at runtime.}
    \item{\code{residuals}:}{A matrix, specifically, the column vector of residuals calculated at runtime.}
    \item{\code{dims}:}{Object of class \code{"character"}.}
    \item{\code{numStats}:}{Numeric; number of observed statistics.}
    \item{\code{dataColumns}:}{Object of class \code{"numeric"}.}
    \item{\code{name}:}{Object of class \code{"character"}.}
    \item{\code{data}:}{Object of class \code{"MxCharOrNumber"}.}
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
Kirkpatrick RM, Pritikin JN, Hunter MD, & Neale MC.  (2021).  Combining structural-equation modeling with genomic-relatedness matrix restricted maximum likelihood in OpenMx.  Behavior Genetics 51: 331-342.  \doi{10.1007/s10519-020-10037-5}

The first software implementation of "GREML":\cr
Yang J, Lee SH, Goddard ME, Visscher PM.  (2011).  GCTA: a tool for genome-wide complex trait analysis.  American Journal of Human Genetics 88: 76-82.  \doi{10.1016/j.ajhg.2010.11.011}

One of the first uses of the acronym "GREML":\cr
Benjamin DJ, Cesarini D, van der Loos MJHM, Dawes CT, Koellinger PD, et al. (2012).  The genetic architecture of economic and political preferences. Proceedings of the National Academy of Sciences 109: 8026-8031. doi: 10.1073/pnas.1120666109

The OpenMx User's guide can be found at \url{https://openmx.ssri.psu.edu/documentation/}.
}

\seealso{
See \code{\link{mxExpectationGREML}()} for creating MxExpectationGREML objects, and for more information generally concerning GREML analyses, including a complete example. More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}
\examples{
showClass("MxExpectationGREML")
}
\keyword{classes}
