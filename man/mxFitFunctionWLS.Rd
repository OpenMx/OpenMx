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

\name{mxFitFunctionWLS}
\alias{mxFitFunctionWLS}
\alias{MxFitFunctionWLS-class}
\alias{print,MxFitFunctionWLS-method}
\alias{show,MxFitFunctionWLS-method}

\title{Create MxFitFunctionWLS Object}

\description{
   This function creates a new MxFitFunctionWLS object.
}

\usage{
mxFitFunctionWLS(type=c('WLS','DWLS','ULS'),
			     allContinuousMethod=c("cumulants", "marginals"),
			     fullWeight=TRUE)
}

\arguments{
   \item{type}{A character string 'WLS' (default), 'DWLS', or 'ULS' for
   weighted, diagonally weighted, or unweighted least squares, respectively}
 \item{allContinuousMethod}{A character string 'cumulants' (default) or
   'marginals'. See Details.}
   \item{fullWeight}{Logical determining if the full weight matrix is
   returned (default). Needed for standard error and quasi-chi-squared
   calculation.}
}

\details{
Fit functions are functions for which free parameter values are optimized such that the value of a cost function is minimized. The mxFitFunctionWLS function computes the weighted least squares difference between the data and the model-implied expectations for the data based on the free parameters and the expectation function (e.g., \link{mxExpectationNormal} or \link{mxExpectationRAM}) selected for the model. 

\strong{Bias and sensitivity to model misspecification}
Both ordinal and continuous data are supported.  A combination of these data  types also succeeds. All three methods ('WLS', 'ULS' and 'DWLS') are unbiased when the model is correct. Full 'WLS', however, is highly sensitive to model misspecification. It can heavily weight the fourth-order moments of the distribution, so small deviations between the observed fourth-order moments and those implied by the model can lead to poor estimates.

\strong{Behavior with all-continuous data}
When only continuous variables are present, the argument \code{allContinuousMethod} dictates how to process the data.

The default, \emph{cumulants} is a good choice for non-normal data. This uses the asymptotically distribution free (ADF) method of Browne (1984) and computes the fourth order \emph{cumulants} for the weight matrix: thus, the name.  The \emph{cumulants} method does not readily handle missing data and does not return weights or summary statistics for the means. It is, however generally fast and ADF up to elliptical distributions. Data computed using cumulants should also be more accurate than via marginals because the whole covariance is a single analytic expression, with no estimation is involved. 

The alternative option, 'marginals', uses methods similar to those used in processing ordinal and joint ordinal-continuous data.  By contrast with cumulats, marginals returns weights and summary statistics for the means.

When data are not all continuous, \code{allContinuousMethod} is ignored.

Usage Notes: 

The results of the optimization can be reported using the \link{summary} function, or accessed directly in the 'output' slot of the resulting model (i.e., modelName$output). Components of the output may be referenced using the \link{Extract} functionality.
}

\value{
Returns a new MxFitFunctionWLS object. One and only one MxFitFunctionWLS object should be included in each model along  with an associated \link{mxExpectationNormal} or \link{mxExpectationRAM} object.
}

\references{
The OpenMx User's guide can be found at http://openmx.ssri.psu.edu/documentation.

Browne, M. W. (1984). Asymptotically distribution-free methods for the analysis of covariance structures. British Journal of Mathematical and Statistical Psychology, 37, p. 62-83.
}

\seealso{
Other fit functions:
\code{\link{mxFitFunctionMultigroup}}, \code{\link{mxFitFunctionML}}, 
\code{\link{mxFitFunctionAlgebra}}, 
\code{\link{mxFitFunctionGREML}}, \code{\link{mxFitFunctionR}}, 
\code{\link{mxFitFunctionRow}}

More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}

\examples{

# Create and fit a model using mxMatrix, mxAlgebra, mxExpectationNormal, and mxFitFunctionWLS

library(OpenMx)

# Simulate some data

x=rnorm(1000, mean=0, sd=1)
y= 0.5*x + rnorm(1000, mean=0, sd=1)
tmpFrame <- data.frame(x, y)
tmpNames <- names(tmpFrame)
wdata <- mxData(tmpFrame, 'raw')

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
