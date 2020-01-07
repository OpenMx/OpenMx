%
%   Copyright 2007-2019 by the individuals mentioned in the source code history
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
mxFitFunctionML(vector = FALSE, rowDiagnostics = FALSE, ..., fellner =
  as.logical(NA), verbose=0L, profileOut=c(),
  rowwiseParallel=as.logical(NA), jointConditionOn = c("auto", "ordinal", "continuous"))
}

\arguments{
   \item{vector}{A logical value indicating whether the objective function result is the likelihood vector.}
   \item{rowDiagnostics}{A logical value indicating whether the row-wise
   results of the objective function should be returned as an attribute
   of the fit function.}
   \item{...}{Not used.  Forces remaining arguments to be specified by name.}
   \item{fellner}{Whether to fully expand the covariance matrix for
     maximum flexibility.}
   \item{verbose}{Level of diagnostic output}
   \item{profileOut}{Character vector naming constant coefficients to
     profile out of the likelihood (sometimes known as REML)}
   \item{rowwiseParallel}{For raw data only, whether to use OpenMP to parallelize the
     evaluation of rows}
   \item{jointConditionOn}{The evaluation strategy when both continuous
   and ordinal data are present.}
}

\details{
Fit functions are functions for which free parameter values are optimized such that the value of a cost function is minimized. The mxFitFunctionML function computes -2*(log likelihood) of the data given the current values of the free parameters and the expectation function (e.g., \link{mxExpectationNormal} or \link{mxExpectationRAM})  selected for the model. 

The 'vector' argument is either TRUE or FALSE, and determines whether the objective function returns a column vector of the likelihoods, or a single -2*(log likelihood) value.

The 'rowDiagnostics' argument is either TRUE or FALSE, and determines whether the row likelihoods are returned as an attribute of the fit function.  Additionally, the squared Mahalanobis distance and the number of observed (non-missing) variables) for each row are returned under the names \code{rowDist} and \code{rowObs}, respectively.  It is sometimes useful to inspect the likelihoods for outliers, diagnostics, or other anomalies.  Each rowwise squared Mahalanobis distance should be chi-squared distributed with degrees of freedom equal to the number of observed variables.  In the case of no missing data, all of the rowwise squared Mahalanobis distances should theoretically be chi-squared distributed with the same degrees of freedom.  In the case of some missing data, the rowwise squared Mahalanobis distances should theoretically be a mixture of chi-squared distributions with mixing proportions equal to the proportions of each number of observed variables.

If there are ordinal data, then only the row likelihoods are returned among the row diagnostics.

When \code{vector=FALSE} and \code{rowDiagnostics=TRUE}, the fit function can be referenced in the model and included in algebras as a scalar. The row likelihoods, row distances, and row observations are then an attribute of the fit function but are not accessible in the model during optimization.  The row likelihoods and other diagnostics are accessible to the user after the model has been run.

By default, \code{jointConditionOn='auto'} and a heuristic will be used
to select the fastest algorithm. Conditioning the continuous data on
ordinal will be superior when there are relatively few unique ordinal
patterns. Otherwise, conditioning the ordinal data on continuous will
perform better when there are relatively many ordinal patterns.

Usage Notes: 

The results of the optimization can be reported using the \link{summary} function, or accessed directly in the 'output' slot of the resulting model (i.e., modelName$output). Components of the output may be referenced using the \link{Extract} functionality.
}

\value{
Returns a new MxFitFunctionML object. One and only one MxFitFunctionML object should be included in each model along  with an associated \link{mxExpectationNormal} or \link{mxExpectationRAM} object.
}

\seealso{
Other fit functions:
\code{\link{mxFitFunctionMultigroup}},  
\code{\link{mxFitFunctionWLS}}, \code{\link{mxFitFunctionAlgebra}}, 
\code{\link{mxFitFunctionGREML}}, \code{\link{mxFitFunctionR}}, 
\code{\link{mxFitFunctionRow}}

More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}

\references{
The OpenMx User's guide can be found at http://openmx.ssri.psu.edu/documentation.
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
