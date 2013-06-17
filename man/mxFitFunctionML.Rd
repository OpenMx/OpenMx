%
%   Copyright 2007-2013 The OpenMx Project
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

\name{mxFitFunctionML}
\alias{mxFitFunctionML}

\title{Create MxFitFunctionML Object}

\description{
   This function creates a new MxFitFunctionML object.
}

\usage{
MxFitFunctionML(vector = FALSE)
}

\arguments{
   \item{vector}{A logical value indicating whether the objective function result is the likelihood vector.}
}

\details{
Fit functions are functions for which free parameter values are optimized such that the value of a cost function is minimized. The mxFitFunctionML function computes -2*(log likelihood) of the data given the current values of the free parameters and the expectation function (e.g., \link{mxExpectationNormal} or \link{mxExpectationRAM})  selected for the model. 

The 'vector' argument is either TRUE or FALSE, and determines whether the objective function returns a column vector of the likelihoods, or a single -2*(log likelihood) value.

Usage Notes: 

The results of the optimization can be reported using the \link{summary} function, or accessed directly in the 'output' slot of the resulting model (i.e., modelName@output). Components of the output may be referenced using the \link{Extract} functionality.
}

\value{
Returns a new MxFitFunctionML object. One and only one MxFitFunctionML object should be included in each model along  with an associated \link{mxExpectationNormal} or \link{mxExpectationRAM} object.
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

expectationFunction <- mxExpectationNormal("expectedCov", "M")

fitFunction <- mxFitFunctionML()

model <- mxModel(A, D, expectedCov, expectationFunction, fitFunction)

}
