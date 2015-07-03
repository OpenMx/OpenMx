%
%   Copyright 2007-2015 The OpenMx Project
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
\alias{MxFitFunctionML-class}
\alias{print,MxFitFunctionML-method}
\alias{show,MxFitFunctionML-method}

\title{Create MxFitFunctionML Object}

\description{
   This function creates a new MxFitFunctionML object.
}

\usage{
mxFitFunctionML(vector = FALSE, rowDiagnostics = FALSE)
}

\arguments{
   \item{vector}{A logical value indicating whether the objective function result is the likelihood vector.}
   \item{rowDiagnostics}{A logical value indicating whether the row-wise results of the objective function should be returned as an attribute of the fit function.}
}

\details{
Fit functions are functions for which free parameter values are optimized such that the value of a cost function is minimized. The mxFitFunctionML function computes -2*(log likelihood) of the data given the current values of the free parameters and the expectation function (e.g., \link{mxExpectationNormal} or \link{mxExpectationRAM})  selected for the model. 

The 'vector' argument is either TRUE or FALSE, and determines whether the objective function returns a column vector of the likelihoods, or a single -2*(log likelihood) value.

The 'rowDiagnostics' arguent is either TRUE or FALSE, and determines whether the row likelihoods are returned as an attribute of the fit function.  It is sometimes useful to inspect the likelihoods for outliers, diagnostics, or other anomalies.

When \code{vector=FALSE} and \code{rowDiagnostics=TRUE}, fitfunction can be referenced in the model and included in algebras as a scalar. The row likelihoods are an attribute of the fit function but are not accessible in the model during optimization.  The row likelihoods are accessible to the user after the model has been run.

Usage Notes: 

The results of the optimization can be reported using the \link{summary} function, or accessed directly in the 'output' slot of the resulting model (i.e., modelName$output). Components of the output may be referenced using the \link{Extract} functionality.
}

\value{
Returns a new MxFitFunctionML object. One and only one MxFitFunctionML object should be included in each model along  with an associated \link{mxExpectationNormal} or \link{mxExpectationRAM} object.
}

\seealso{\link{mxFitFunctionMultigroup} for multiple group models and \link{mxFitFunctionAlgebra} for user-defined fit functions.}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{

# Create and fit a model using mxMatrix, mxAlgebra, mxExpectationNormal, and mxFitFunctionML

library(OpenMx)

# Simulate some data

x=rnorm(1000, mean=0, sd=1)
y= 0.5*x + rnorm(1000, mean=0, sd=1)
tmpFrame <- data.frame(x, y)
tmpNames <- names(tmpFrame)

# Define the matrices

M <- mxMatrix(type = "Full", nrow = 1, ncol = 2, values=c(0,0), 
              free=c(TRUE,TRUE), labels=c("Mx", "My"), name = "M")
S <- mxMatrix(type = "Full", nrow = 2, ncol = 2, values=c(1,0,0,1), 
              free=c(TRUE,FALSE,FALSE,TRUE), labels=c("Vx", NA, NA, "Vy"), name = "S")
A <- mxMatrix(type = "Full", nrow = 2, ncol = 2, values=c(0,1,0,0), 
              free=c(FALSE,TRUE,FALSE,FALSE), labels=c(NA, "b", NA, NA), name = "A")
I <- mxMatrix(type="Iden", nrow=2, ncol=2, name="I")

# Define the expectation

expCov <- mxAlgebra(solve(I-A) \%*\% S \%*\% t(solve(I-A)), name="expCov")
expFunction <- mxExpectationNormal(covariance="expCov", means="M", dimnames=tmpNames)

# Choose a fit function

fitFunction <- mxFitFunctionML(rowDiagnostics=TRUE)
# also return row likelihoods, even though the fit function
#  value is still 1x1

# Define the model

tmpModel <- mxModel(model="exampleModel", M, S, A, I, expCov, expFunction, fitFunction, 
                    mxData(observed=tmpFrame, type="raw"))

# Fit the model and print a summary

tmpModelOut <- mxRun(tmpModel)
summary(tmpModelOut)

fitResOnly <- mxEval(fitfunction, tmpModelOut)
attributes(fitResOnly) <- NULL
fitResOnly

# Look at the row likelihoods alone
fitLikeOnly <- attr(mxEval(fitfunction, tmpModelOut), 'likelihoods')
head(fitLikeOnly)

}
