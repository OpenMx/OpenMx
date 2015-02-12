\name{mxFitFunctionGREML}
\alias{mxFitFunctionGREML}

\title{
Create MxFitFunctionGREML Object
}
\description{
This function creates a new \code{\link{MxFitFunctionGREML}} object.
}
\usage{
mxFitFunctionGREML(casesToDrop=integer(0), dropNAfromV=TRUE)
}

\arguments{
  \item{casesToDrop}{Indices for the rows and columns to be internally dropped from the 'V' matrix (and its derivatives, if applicable), usually because those same rows have been dropped from 'X' and 'y' due to missing data.  Defaults to an integer vector of length zero.}
  \item{dropNAfromV}{Logical; should the rows and columns indexed by argument \code{casesToDrop} actually be dropped at runtime?  Defaults to \code{TRUE}.}
}
\details{
Rows and columns are internally dropped only if argument \code{casesToDrop} has nonzero length and argument \code{dropNAfromV} is \code{TRUE}.  The point of these two arguments is that, frequently, the 'V' matrix will be an \code{\link{MxAlgebra}}, and therefore not be computed until runtime.  In contrast, 'y' and 'X' will usually be \code{\link{MxMatrix}} objects and can easily have rows with missing data trimmed out of them prior to runtime. 
}
\value{
Returns a new object of class \code{\link{MxFitFunctionGREML}}.
}
\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\seealso{
See \code{\link{mxGREMLStarter}()} for a related helper function.  See \code{\link{MxFitFunctionGREML}} for the S4 class created by \code{mxFitFunctionGREML()}.  For more information generally concerning GREML analyses, see \code{\link{mxExpectationGREML}()}.  More information about the OpenMx package may be found \link[=OpenMx]{here}.
}
\examples{
gff <- mxFitFunctionGREML()
}
