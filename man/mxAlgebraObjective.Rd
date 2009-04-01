\name{mxAlgebraObjective}
\alias{mxAlgebraObjective}

\title{Function to Create MxAlgebraObjective Object}

\description{
   This function creates a new MxAlgebraObjective object.
}

\usage{
mxAlgebraObjective(algebra, name = NA)
}

\arguments{
   \item{algebra}{A character string indicating the name of an \link{MxAlgebra} or \link{MxMatrix} object to use for optimization.}
   \item{name}{An optional character string indicating the name of this object.}
}

\details{
Objective functions are functions for which free parameter values are chosen such that the value of the objective function is minimized. While the other objective functions in OpenMx are packaged with a function to be optimized (i.e., maximum likelihood), the mxAlgebraObjective function uses the referenced \link{MxAlgebra} or \link{MxMatrix} object as the function to be minimized.

If no free parameters are included in any part of the 'algebra' argument, the mxAlgebraObjective works as an algebra evaluator. If the 'algebra' argument is given an \link{MxMatrix} object with no free parameters, the original object is returned. If the 'algebra' argument is given an \link{MxAlgebra} object whose referenced matrices and algebras have no free parameters, the mxAlgebraObjective function carries out the operations defined in the 'expression' argument of that \link{MxAlgebra} object.

If free parameters are included in the objects references in the 'algebra' argument, the mxAlgebraObjective minimizes the algebra defined in the 'algebra' argument. More than one free parameter may be included in the algebra, but the mxAlgebraObjective should be specified such that the objective function evaluates to a 1x1 matrix or scalar value. If the contents of the 'algebra' argument evaluates to a non-scalar matrix, only the first entry of that matrix is minimized.

To evaluate, place MxAlgebraObjective objects, referenced \link{MxAlgebra} and \link{MxMatrix} objects, and optional \link{MxBounds} and \link{MxConstraint} objects in an \link{MxModel} object. This model may then be evaluated using the \link{mxRun} function. The results of the optimization can be found in the 'output' slot of the resulting model, and may be referenced using the \link{Extract} functionality.
}

\value{
Returns a new MxAlgebraObjective object. MxAlgebraObjective objects should be included with models with referenced \link{MxAlgebra} and \link{MxMatrix} objects.
}

\references{
The Mx manual can be found at http://www.vcu.edu/mx/documentation.html.
}

\examples{
#Create a matrix 'A' with no free parameters
A <- mxMatrix('Full', nrow = 3, ncol = 3, values = 1:9, name = 'A')

#Create an algebra 'B', which defines the expression A+A
B <- mxAlgebra(A + A, name = 'B')

#Define the objective function for algebra 'B'
objective <- mxAlgebraObjective('B')

#Place the algebra, its associated matrix and its objective function in a model
model <- mxModel(A, B, objective)

#Evalulate the algebra
modelRun <- mxRun(model)

#View the results
modelRun@output
}