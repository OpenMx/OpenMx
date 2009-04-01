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
   \item{means}{A character string indicating the name of expected means algebra. Optional.}
   \item{name}{A character string indicating the
   		name of this object. Optional.}
}

\value{
	Returns a new MxMLObjective object.
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
