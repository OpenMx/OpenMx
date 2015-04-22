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

\name{mxDataWLS}
\alias{mxDataWLS}

\title{Create MxData Object for Least Squares (WLS, DLS, ULS) Analyses}

\description{
   This function creates a new \link{MxData} object of type 
   \dQuote{ULS} (unweighted least squares), \dQuote{WLS} (weighted least squares) 
   or \dQuote{DLS} (diagonally-weighted least squares). The appropriate
   fit function to include with these models is \code{\link{mxFitFunctionWLS}}
}

\usage{
   mxDataWLS(data, type = "WLS", useMinusTwo = TRUE, returnInverted = TRUE, 
    debug = FALSE, fullWeight = TRUE)
}

\arguments{
   \item{data}{A matrix or data.frame which provides raw data to be used for WLS.}
   \item{type}{A character string 'WLS' (default), 'DLS', or 'ULS' for weighted, diagonal, or unweighted least squares}
   \item{useMinusTwo}{Logical indicating whether to use -2LL (default) or -LL.}
   \item{returnInverted}{Logical indicating whether to return the information matrix (default) or the covariance matrix.}
   \item{debug}{Logical to set debugging on or off (default)}
   \item{fullWeight}{Logical determining if the full weight matrix is returned (default). Needed for standard error and quasi-chi-squared calculation.}
}

\details{
The mxDataWLS function creates an \link{MxData} object, which can be used in
\link{MxModel} objects.  This function takes raw data and returns 
an \code{MxData} object to be used in a model to fit with weighted least squares.

Both Ordinal and continuous data are supported.  A combination of these data 
types succeeds without error, but when using 'WLS' or 'DLS' the answers
appear incorrect.  The 'ULS' estimates for joint ordinal and continuous data
appear accurate.  Consequently, do not use this function for joint 
problems unless \code{type='ULS'}.
}

\value{
    Returns a new \link{MxData} object.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\seealso{
\link{mxFitFunctionWLS}.  \link{MxData} for the S4 class created by mxData. \link{matrix} and \link{data.frame} for objects which may be entered as arguments in the \sQuote{observed} slot. More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}

\examples{

# Create and fit a model using mxMatrix, mxAlgebra, mxExpectationNormal, and mxFitFunctionWLS

library(OpenMx)

# Simulate some data

x=rnorm(1000, mean=0, sd=1)
y= 0.5*x + rnorm(1000, mean=0, sd=1)
tmpFrame <- data.frame(x, y)
tmpNames <- names(tmpFrame)
wdata <- mxDataWLS(tmpFrame)

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

fitFunction <- mxFitFunctionWLS()

# Define the model

tmpModel <- mxModel(model="exampleModel", S, A, I, expCov, expFunction, fitFunction, 
                    wdata)

# Fit the model and print a summary

tmpModelOut <- mxRun(tmpModel)
summary(tmpModelOut)
}
