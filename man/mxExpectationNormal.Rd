%
%   Copyright 2007-2021 by the individuals mentioned in the source code history
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

\name{mxExpectationNormal}
\alias{mxExpectationNormal}
\alias{MxExpectationNormal-class}
\alias{show,MxExpectationNormal-method}
\alias{print,MxExpectationNormal-method}

\title{Create MxExpectationNormal Object}

\description{
   This function creates an MxExpectationNormal object.
}

\usage{
mxExpectationNormal(covariance, means=NA, dimnames = NA, thresholds = NA,
                    threshnames = dimnames, ...,
  discrete = as.character(NA), discreteSpec=NULL)
}

\arguments{
   \item{covariance}{A character string indicating the name of the expected covariance algebra.}
   \item{means}{A character string indicating the name of the expected means algebra.}
   \item{dimnames}{An optional character vector to be assigned to the dimnames of the covariance and means algebras.}
   \item{thresholds}{An optional character string indicating the name of the thresholds matrix.}
   \item{threshnames}{An optional character vector to be assigned to the column names of the thresholds matrix.}
   \item{...}{Not used.  Forces remaining arguments to be specified by name.}
   \item{discrete}{An optional character string indicating the name of the discrete matrix.}
   \item{discreteSpec}{An optional matrix containing maximum counts and
   model IDs.}
}

\details{
Expectation functions define the way that model expectations are calculated. The mxExpectationNormal function uses the algebra defined by the 'covariance' and 'means' arguments to define the expected covariance and means under the assumption of multivariate normality. The 'covariance' argument takes an \link{MxAlgebra} object, which defines the expected covariance of an associated \link{MxData} object. The 'means' argument takes an \link{MxAlgebra} object, which defines the expected means of an associated \link{MxData} object. The 'dimnames' arguments takes an optional character vector.  If this argument is not a single NA, then this vector is used to assign the dimnames of the means vector as well as the row and columns dimnames of the covariance matrix.

thresholds: The name of the thresholds matrix. When needed (for modelling ordinal data), this matrix should be created using  \link{mxMatrix}(). The thresholds matrix must have as many columns as there are ordinal variables in the model, and number of rows equal to one fewer than the maximum number of levels found in the ordinal variables. The starting values  of this matrix must also be set to reasonable values. Fill each column with a set of ordered start thresholds, one for each level of this column's factor levels minus 1. These thresholds may be free if you wish them to be estimated, or fixed. The unused rows in each column, if any, can be set to any value including NA.

threshnames: A character vector consisting of the variables in the thresholds matrix, i.e., the names of ordinal variables in a model. This is necessary for OpenMx to map the thresholds matrix columns onto the variables in your data. If you set the \link{dimnames} of the columns in the thresholds matrix then threshnames is not needed.

Usage Notes:
dimnames must be supplied where the matrices referenced by the covariance and means algebras are not themselves labeled. Failure to do so leads to an error noting that the covariance or means matrix associated with the FIML objective does not contain dimnames.

mxExpectationNormal evaluates with respect to an \link{MxData} object. The \link{MxData} object need not be referenced in the mxExpectationNormal function, but must be included in the \link{MxModel} object. When the 'type' argument in the associated \link{MxData} object is equal to 'raw', missing values are permitted in the associated \link{MxData} object.

To evaluate, place an mxExpectationNormal object, the \link{mxData} object for which the expected covariance approximates, referenced \link{MxAlgebra} and \link{MxMatrix} objects, optional \link{MxBounds} or \link{MxConstraint} objects, and an mxFitFunction such as \link{mxFitFunctionML} in an \link{MxModel} object. This model may then be evaluated using the \link{mxRun} function.

The results of the optimization can be reported using the \link{summary} function, or accessed directly in the 'output' slot of the resulting model (i.e., modelName$output). Components of the output may be referenced using the \link{Extract} functionality.
}

\value{
Returns an MxExpectationNormal object.
}

\references{
The OpenMx User's guide can be found at \url{https://openmx.ssri.psu.edu/documentation/}.
}

\examples{

# Create and fit a model using mxMatrix, mxAlgebra,
#  mxExpectationNormal, and mxFitFunctionML

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
              free=c(TRUE,FALSE,FALSE,TRUE), labels=c("Vx", NA, NA, "Vy"),
              name = "S")
A <- mxMatrix(type = "Full", nrow = 2, ncol = 2, values=c(0,1,0,0),
              free=c(FALSE,TRUE,FALSE,FALSE), labels=c(NA, "b", NA, NA),
              name = "A")
I <- mxMatrix(type="Iden", nrow=2, ncol=2, name="I")

# Define the expectation

expCov <- mxAlgebra(solve(I-A) \%*\% S \%*\% t(solve(I-A)), name="expCov")
expFunction <- mxExpectationNormal(covariance="expCov", means="M",
	dimnames=tmpNames)

# Choose a fit function

fitFunction <- mxFitFunctionML()

# Define the model

tmpModel <- mxModel(model="exampleModel", M, S, A, I,
                    expCov, expFunction, fitFunction,
                    mxData(observed=tmpFrame, type="raw"))

# Fit the model and print a summary

tmpModelOut <- mxRun(tmpModel)
summary(tmpModelOut)


}
