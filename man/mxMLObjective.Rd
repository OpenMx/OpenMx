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

\name{mxMLObjective}
\alias{mxMLObjective}

\title{DEPRECATED: Create MxMLObjective Object}

\description{
WARNING: Objective functions have been deprecated as of OpenMx 2.0.  
   
Please use mxExpectationNormal() and mxFitFunctionML() instead.  As a temporary workaround, mxMLObjective returns a list containing an MxExpectationNormal object and an MxFitFunctionML object.
}

\usage{
mxMLObjective(covariance, means = NA, dimnames = NA, thresholds = NA)
All occurrences of

mxMLObjective(covariance, means = NA, dimnames = NA, thresholds = NA)

Should be changed to

mxExpectationNormal(covariance, means = NA, dimnames = NA, thresholds = NA, threshnames = dimnames)
mxFitFunctionML(vector = FALSE)
}


\arguments{
   \item{covariance}{A character string indicating the name of the expected covariance algebra.}
   \item{means}{An optional character string indicating the name of the expected means algebra.}
   \item{dimnames}{An optional character vector to be assigned to the dimnames of the covariance and means algebras.}   
   \item{thresholds}{An optional character string indicating the name of the thresholds matrix.}
}

\details{
NOTE: THIS DESCRIPTION IS DEPRECATED.  Please change to using \link{mxExpectationNormal} and \link{mxFitFunctionML} as shown in the example below.

Objective functions are functions for which free parameter values are chosen such that the value of the objective function is minimized. The mxMLObjective function uses full-information maximum likelihood to provide maximum likelihood estimates of free parameters in the algebra defined by the 'covariance' argument given the covariance of an \link{MxData} object. The 'covariance' argument takes an \link{MxAlgebra} object, which defines the expected covariance of an associated \link{MxData} object. The 'dimnames' arguments takes an optional character vector.  If this argument is not a single NA, then this vector be assigned to be the dimnames of the means vector, and the row and columns dimnames of the covariance matrix.

mxMLObjective evaluates with respect to an \link{MxData} object. The \link{MxData} object need not be referenced in the mxMLObjective function, but must be included in the \link{MxModel} object. mxMLObjective requires that the 'type' argument in the associated \link{MxData} object be equal to 'cov', 'cov', or 'sscp'. The 'covariance' argument of this function evaluates with respect to the 'matrix' argument of the associated \link{MxData} object, while the 'means' argument of this function evaluates with respect to the 'vector' argument of the associated \link{MxData} object. The 'means' and 'vector' arguments are optional in both functions. If the 'means' argument is not specified (NA), the optional 'vector' argument of the \link{MxData} object is ignored. If the 'means' argument is specified, the associated \link{MxData} object should specify a 'means' argument of equivalent dimension as the 'means' algebra.

dimnames must be supplied where the matrices referenced by the covariance and means algebras are not themselves labeled. Failure to do so leads to an error noting that the covariance or means matrix associated with the ML objective does not contain dimnames.

To evaluate, place MxMLObjective objects, the \link{mxData} object for which the expected covariance approximates, referenced \link{MxAlgebra} and \link{MxMatrix} objects, and optional \link{MxBounds} and \link{MxConstraint} objects in an \link{MxModel} object. This model may then be evaluated using the \link{mxRun} function. The results of the optimization can be found in the 'output' slot of the resulting model, or using the \link{mxEval} function.
}

\value{
Returns a list containing an MxExpectationNormal object and an MxFitFunctionML object. 
}

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

S <- mxMatrix(type = "Full", nrow = 2, ncol = 2, values=c(1,0,0,1), 
              free=c(TRUE,FALSE,FALSE,TRUE), labels=c("Vx", NA, NA, "Vy"), name = "S")
A <- mxMatrix(type = "Full", nrow = 2, ncol = 2, values=c(0,1,0,0), 
              free=c(FALSE,TRUE,FALSE,FALSE), labels=c(NA, "b", NA, NA), name = "A")
I <- mxMatrix(type="Iden", nrow=2, ncol=2, name="I")

# Define the expectation

expCov <- mxAlgebra(solve(I-A) \%*\% S \%*\% t(solve(I-A)), name="expCov")
expFunction <- mxExpectationNormal(covariance="expCov", dimnames=tmpNames)

# Choose a fit function

fitFunction <- mxFitFunctionML()

# Define the model

tmpModel <- mxModel("exampleModel", S, A, I, expCov, expFunction, fitFunction, 
                    mxData(observed=cov(tmpFrame), type="sscp", numObs=dim(tmpFrame)[1]))

# Fit the model and print a summary

tmpModelOut <- mxRun(tmpModel)
summary(tmpModelOut)


}
