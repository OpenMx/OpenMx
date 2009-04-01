\name{mxMLObjective}
\alias{mxMLObjective}

\title{Function To Create MxMLObjective Object}

\description{
   This function creates a new MxMLObjective object.
}

\usage{
mxMLObjective(covariance, means = NA, name = NA)
}


\arguments{
   \item{covariance}{A character string indicating the name of the expected covariance algebra.}
   \item{means}{An optional character string indicating the name of expected means algebra.}
   \item{name}{An optional character string indicating the
   		name of this object.}
}

\details{
Objective functions are functions for which free parameter values are chosen such that the value of the objective function is minimized. The mxMLObjective function uses full-information maximum likelihood to provide maximum likelihood estimates of free parameters in the algebra defined by the 'covariance' argument given the covariance of an \link{MxData} object. The 'covariance' argument takes an \link{MxAlgebra} object, which defines the expected covariance of an associated \link{MxData} object.

mxMLObjective evaluates with respect to an \link{MxData} object. The \link{MxData} object need not be referenced in the mxMLObjective function, but must be included in the \link{MxModel} object. mxMLObjective requires that the 'type' argument in the associated \link{MxData} object be equal to 'cov', 'cov', or 'sscp'. The 'covariance' argument of this function evaluates with respect to the 'matrix' argument of the associated \link{MxData} object, while the 'means' argument of this function evaluates with respect to the 'vector' argument of the associated \link{MxData} object. The 'means' and 'vector' arguments are optional in both functions. If the 'means' argument is not specified (NA), the optional 'vector' argument of the \link{MxData} object is ignored. If the 'means' argument is specified, the associated \link{MxData} object should specify a 'means' argument of equivalent dimension as the 'means' algebra.

To evaluate, place MxMLObjective objects, the \link{mxData} object for which the expected covariance approximates, referenced \link{MxAlgebra} and \link{MxMatrix} objects, and optional \link{MxBounds} and \link{MxConstraint} objects in an \link{MxModel} object. This model may then be evaluated using the \link{mxRun} function. The results of the optimization can be found in the 'output' slot of the resulting model, and may be referenced using the \link{Extract} functionality.
}

\value{
	Returns a new MxMLObjective object. MxMLObjective objects should be included with models with referenced \link{MxAlgebra}, \link{MxData} and \link{MxMatrix} objects.
}

\references{
The Mx manual can be found at http://www.vcu.edu/mx/documentation.html.
}

\examples{

A <- mxMatrix(values = 0.5, nrow = 2, ncol = 1, 
	free = TRUE, name = "A")

D <- mxMatrix(type = "Diag", values = c(0, 0.5), 
	specification = c(0, NA), nrow = 2, free = TRUE, 
	name = "D")

expectedCov <- mxAlgebra(A \%*\% t(A) + D, "expectedCov")

observedCov <- mxData(matrix(c(1.2, 0.8, 0.8, 1.3)), 
	'cov', numObs = 150, name = 'observedCov')

objective <- mxMLObjective(covariance = "expectedCov",
	name = "objective")

model <- mxModel(A, D, expectedCov, objective, observedCov)
}
