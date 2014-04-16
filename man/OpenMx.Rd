%
%   Copyright 2007-2014 The OpenMx Project
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
OpenMx solves algebra optimization and statistical estimation problems using matrix algebra. The OpenMx library defines a set of S4 classes and functions used to create them. The majority of these classes are used as arguments in models, which may include data, matrices, algebras, bounds and constraints. Models most often have an expectation function (e.g., mxExpectationNormal) to calculate the expectations for the model.  Each model has one and only one fit function which may either be built-in (e.g., \link{mxFitFunctionML} ) or user-defined fit functions (e.g., \link{mxFitFunctionR} or \link{mxFitFunctionAlgebra}). The resulting model can then be optimized using \link{mxRun}, resulting in parameter estimates, algebra evaluation, and output for additional models.

The main OpenMx functions are: \link{mxAlgebra}, \link{mxBounds}, \link{mxCI}, \link{mxConstraint}, \link{mxData}, \link{mxMatrix}, \link{mxModel}, and \link{mxPath}; the expectation functions such as \link{mxExpectationNormal}, \link{mxExpectationRAM}, \link{mxExpectationLISREL}, and \link{mxExpectationStateSpace}; and the fit functions such as \link{mxFitFunctionML}, \link{mxFitFunctionAlgebra}, \link{mxFitFunctionRow} and \link{mxFitFunctionR}. These latter functions create the objectives for model estimation. Models are run (estimated) using \link{mxRun}.

These functions create the following Objects: \link{MxModel}s, which can contain \link[=MxPath-class]{MxPath}s, \link{MxMatrix} objects, and \link{MxData} as well as \link{MxAlgebra}, \link{MxBounds}, \link{MxCI}s (confidence intervals) and \link{MxConstraint}s.

OpenMx comes with several useful datasets built-in. Access them using data(package="OpenMx")

The main learning material is online, but functions \link{mxRun}, \link{mxModel}, \link{mxMatrix} all have runnable examples to get you started.

}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}
