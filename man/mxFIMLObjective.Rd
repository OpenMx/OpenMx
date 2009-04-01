\name{mxFIMLObjective}
\alias{mxFIMLObjective}

\title{Function To Create MxFIMLObjective Object}

\description{
   This function creates a new MxFIMLObjective object.
}

\usage{
mxFIMLObjective(covariance, means = NA, name = NA)
}


\arguments{
   \item{covariance}{A character string indicating the name of the expected covariance algebra.}
   \item{means}{An optional character string indicating the name of expected means algebra.}
   \item{name}{An optional character string indicating the
   		name of this object.}
}

\details{
Objective functions are functions for which free parameter values are chosen such that the value of the objective function is minimized. The mxFIMLObjective function uses full-information maximum likelihood to provide maximum likelihood estimates of free parameters in the algebra defined by the 'covariance' argument given the covariance of an \link{MxData} object. 

To evaluate, place MxFIMLObjective objects, the \link{mxData} object for which the expected covariance approximates, referenced \link{MxAlgebra} and \link{MxMatrix} objects, and optional \link{MxBounds} and \link{MxConstraint} objects in an \link{MxModel} object. The \link{MxData} object need not be referenced in the mxFIMLObjective function, but must be included in the \link{MxModel} object. This model may then be evaluated using the \link{mxRun} function. The results of the optimization can be found in the 'output' slot of the resulting model, and may be referenced using the \link{Extract} functionality.
}

\value{
Returns a new MxFIMLObjective object. MxFIMLObjective objects should be included with models with referenced \link{MxAlgebra}, \link{MxData} and \link{MxMatrix} objects.
}

\examples{

A <- mxMatrix(values = 0.5, nrow = 2, ncol = 1, 
	free = TRUE, name = "A")

D <- mxMatrix(type = "Diag", values = c(0, 0.5), 
	specification = c(0, NA), nrow = 2, free = TRUE, 
	name = "D")

expectedCov <- mxAlgebra(A \%*\% t(A) + D, "expectedCov")

objective <- mxFIMLObjective("expectedCov")

model <- mxModel(A, D, expectedCov, objective)

}