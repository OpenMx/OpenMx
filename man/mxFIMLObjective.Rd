\name{mxFIMLObjective}
\alias{mxFIMLObjective}

\title{Function To Create MxFIMLObjective Object}

\description{
   This function creates a new MxFIMLObjective object.
}

\usage{
mxFIMLObjective(covariance, means, dimnames = NA, thresholds = NA)
}


\arguments{
   \item{covariance}{A character string indicating the name of the expected covariance algebra.}
   \item{means}{A character string indicating the name of the expected means algebra.}
   \item{dimnames}{An optional character vector to be assigned to the dimnames of the covariance and means algebras.}
   \item{thresholds}{An optional character string indicating the name of the thresholds matrix.}
}

\details{
Objective functions are functions for which free parameter values are chosen such that the value of the objective function is minimized. The mxFIMLObjective function uses full-information maximum likelihood to provide maximum likelihood estimates of free parameters in the algebra defined by the 'covariance' and 'means' arguments. The 'covariance' argument takes an \link{MxAlgebra} object, which defines the expected covariance of an associated \link{MxData} object. The 'means' argument takes an \link{MxAlgebra} object, which defines the expected means of an associated \link{MxData} object. The 'dimnames' arguments takes an optional character vector.  If this argument is not a single NA, then this vector be assigned to be the dimnames of the means vector, and the row and columns dimnames of the covariance matrix.

mxFIMLObjective evaluates with respect to an \link{MxData} object. The \link{MxData} object need not be referenced in the mxFIMLObjective function, but must be included in the \link{MxModel} object. mxFIMLObjective requires that the 'type' argument in the associated \link{MxData} object be equal to 'raw'. Missing values are permitted in the associated \link{MxData} object.

To evaluate, place MxFIMLObjective objects, the \link{mxData} object for which the expected covariance approximates, referenced \link{MxAlgebra} and \link{MxMatrix} objects, and optional \link{MxBounds} and \link{MxConstraint} objects in an \link{MxModel} object. This model may then be evaluated using the \link{mxRun} function. The results of the optimization can be found in the 'output' slot of the resulting model, and may be referenced using the \link{Extract} functionality.
}

\value{
Returns a new MxFIMLObjective object. MxFIMLObjective objects should be included with models with referenced \link{MxAlgebra}, \link{MxData} and \link{MxMatrix} objects.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{

A <- mxMatrix(values = 0.5, nrow = 2, ncol = 1, 
	free = TRUE, name = "A")

D <- mxMatrix(type = "Diag", values = c(0, 0.5), 
	free = c(FALSE, TRUE), nrow = 2, name = "D")
	
M <- mxMatrix(type = "Zero", nrow = 1, ncol = 2, name = "M")

expectedCov <- mxAlgebra(A \%*\% t(A) + D, "expectedCov")

objective <- mxFIMLObjective("expectedCov", "M")

model <- mxModel(A, D, expectedCov, objective)

}
