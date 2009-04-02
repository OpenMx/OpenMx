\name{OpenMx}
\alias{OpenMx}

\title{OpenMx: Package for Matrix Algebra Optimization}

\description{OpenMx is a package for structural equation modeling, matrix algebra optimization and other statistical estimation problems. 
}

\details{
OpenMx is a package for algebra optimization and statistical estimation problems using matrix algebra. The OpenMx library defines a set of S4 classes and functions used to create them. The majority of these classes are used as arguments in models, which may include data, matrices, algebras, bounds and constraints. These models are then paired with objective functions, either existing (maximum likelihood, FIML) or user-defined with included algebra functions. These models can then be optimized, resulting in parameter estimation, algebra evaluation, and output for additional models.

Objects used or created by OpenMx may be of the following classes: \link{MxAlgebra}, \link{MxBounds}, \link{MxConstraint}, \link{MxData}, \link{MxMatrix}, \link{MxModel}, and \link{MxPath}. Objects of these classes may be created by the following OpenMx functions: \link{mxAlgebra}, \link{mxBounds}, \link{mxConstraint}, \link{mxData}, \link{mxMatrix}, \link{mxModel}, and \link{mxPath}. The functions \link{mxAlgebraObjective}, \link{mxFIMLObjective}, \link{mxMLObjective} and \link{mxRAMObjective} create objective functions for model estimation. Models which include objective functions may be estimated using the \link{mxRun} function.
}

\references{
The Mx manual may be found at http://www.vcu.edu/mx/documentation.html.
}