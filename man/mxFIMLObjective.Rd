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
   \item{means}{A character string indicating the name of expected means algebra. Optional.}
   \item{name}{A character string indicating the
   		name of this object. Optional.}
}

\value{
	Returns a new MxFIMLObjective object.
}

\examples{

A <- mxMatrix(values = 0.5, nrow = 2, ncol = 1, 
	free = TRUE, name = "A")

D <- mxMatrix(type = "Diag", values = c(0, 0.5), 
	specification = c(0, NA), nrow = 2, free = TRUE, 
	name = "D")

expectedCov <- mxAlgebra(A %*% t(A) + D, "expectedCov")

objective <- mxFIMLObjective("expectedCov")

model <- mxModel(A, D, expectedCov, objective)

}