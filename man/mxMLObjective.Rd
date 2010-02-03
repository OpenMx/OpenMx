\name{mxMLObjective}
\alias{mxMLObjective}

\title{Create MxMLObjective Object}

\description{
   This function creates a new MxMLObjective object.
}

\usage{
mxMLObjective(covariance, means = NA, dimnames = NA, thresholds = NA)
}


\arguments{
   \item{covariance}{A character string indicating the name of the expected covariance algebra.}
   \item{means}{An optional character string indicating the name of the expected means algebra.}
   \item{dimnames}{An optional character vector to be assigned to the dimnames of the covariance and means algebras.}   
   \item{thresholds}{An optional character string indicating the name of the thresholds matrix.}
}

\details{
Objective functions are functions for which free parameter values are chosen such that the value of the objective function is minimized. The mxMLObjective function uses full-information maximum likelihood to provide maximum likelihood estimates of free parameters in the algebra defined by the 'covariance' argument given the covariance of an \link{MxData} object. The 'covariance' argument takes an \link{MxAlgebra} object, which defines the expected covariance of an associated \link{MxData} object. The 'dimnames' arguments takes an optional character vector.  If this argument is not a single NA, then this vector be assigned to be the dimnames of the means vector, and the row and columns dimnames of the covariance matrix.

mxMLObjective evaluates with respect to an \link{MxData} object. The \link{MxData} object need not be referenced in the mxMLObjective function, but must be included in the \link{MxModel} object. mxMLObjective requires that the 'type' argument in the associated \link{MxData} object be equal to 'cov', 'cov', or 'sscp'. The 'covariance' argument of this function evaluates with respect to the 'matrix' argument of the associated \link{MxData} object, while the 'means' argument of this function evaluates with respect to the 'vector' argument of the associated \link{MxData} object. The 'means' and 'vector' arguments are optional in both functions. If the 'means' argument is not specified (NA), the optional 'vector' argument of the \link{MxData} object is ignored. If the 'means' argument is specified, the associated \link{MxData} object should specify a 'means' argument of equivalent dimension as the 'means' algebra.

dimnames must be supplied where the matrices referenced by the covariance and means algebras are not themselves labeled. Failure to do so leads to an error noting that the covariance or means matrix associated with the ML objective does not contain dimnames.

To evaluate, place MxMLObjective objects, the \link{mxData} object for which the expected covariance approximates, referenced \link{MxAlgebra} and \link{MxMatrix} objects, and optional \link{MxBounds} and \link{MxConstraint} objects in an \link{MxModel} object. This model may then be evaluated using the \link{mxRun} function. The results of the optimization can be found in the 'output' slot of the resulting model, or using the \link{mxEval} function.
}

\value{
	Returns a new MxMLObjective object. MxMLObjective objects should be included with models with referenced \link{MxAlgebra}, \link{MxData} and \link{MxMatrix} objects.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{

vars <- c('x','y')

A <- mxMatrix(values = 0.5, nrow = 2, ncol = 1,
   free = TRUE, name = "A")
D <- mxMatrix(type = "Diag", values = c(0, 0.5),
   free = c(FALSE, TRUE), nrow = 2, name = "D")

expectedCov <- mxAlgebra(A \%*\% t(A) + D, "expectedCov")
observedCov <- mxData(matrix(c(1.2, 0.8, 0.8, 1.3),
   nrow = 2, ncol = 2, dimnames = list(vars,vars)), 'cov', numObs = 150)

objective <- mxMLObjective(covariance = "expectedCov", dimnames = vars)

model <- mxModel("mxMLObjective example", A, D, expectedCov, objective, observedCov)

\dontrun{summary(mxRun(model))}

}
