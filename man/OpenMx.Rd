%
%   Copyright 2007-2010 The OpenMx Project
%
%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License.
%   You may obtain a copy of the License at
% 
%        http://www.apache.org/licenses/LICENSE-2.0
% 
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.

\name{OpenMx}
\alias{OpenMx}

\title{OpenMx: Package for Matrix Algebra Optimization}

\description{OpenMx is a package for structural equation modeling, matrix algebra optimization and other statistical estimation problems. 
}

\details{
OpenMx is a package for algebra optimization and statistical estimation problems using matrix algebra. The OpenMx library defines a set of S4 classes and functions used to create them. The majority of these classes are used as arguments in models, which may include data, matrices, algebras, bounds and constraints. These models are then paired with objective functions, either existing (maximum likelihood, FIML) or user-defined with included algebra functions. These models can then be optimized, resulting in parameter estimation, algebra evaluation, and output for additional models.

Objects used or created by OpenMx may be of the following classes: \link{MxAlgebra}, \link{MxBounds}, \link{MxCI}, \link{MxConstraint}, \link{MxData}, \link{MxMatrix}, \link{MxModel}, and \link{MxPath}. Objects of these classes may be created by the following OpenMx functions: \link{mxAlgebra}, \link{mxBounds}, \link{mxCI}, \link{mxConstraint}, \link{mxData}, \link{mxMatrix}, \link{mxModel}, and \link{mxPath}. The functions \link{mxAlgebraObjective}, \link{mxFIMLObjective}, \link{mxMLObjective} and \link{mxRAMObjective} create objective functions for model estimation. Models which include objective functions may be estimated using the \link{mxRun} function.

OpenMx comes with several useful datasets built-in. Access them using data(package="OpenMx")
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}
