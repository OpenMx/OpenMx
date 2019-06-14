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
   \item{type}{A character string 'WLS' (default), 'DWLS', or 'ULS' for weighted, diagonally weighted, or unweighted least squares, respectively}
   \item{allContinuousMethod}{A character string 'cumulants' (default) or 'marginals'. See Details.}
   \item{fullWeight}{Logical determining if the full weight matrix is returned (default). Needed for standard error and quasi-chi-squared calculation.}
}

\details{
As with other fit functions, \code{mxFitFunctionWLS} optimizes free parameter values such that the value of a cost function is minimized. For \code{mxFitFunctionWLS}, this cost function is the weighted least squares difference between the data and the model-implied expectations for the data based on the free parameters and the expectation function (e.g., \code{\link{mxExpectationNormal}} or \code{\link{mxExpectationRAM}}) selected for the model. 

\strong{Bias and sensitivity to model misspecification}
Both ordinal and continuous data are supported, as well as combinations of these data types. All three methods ('WLS', 'ULS' and 'DWLS') are unbiased when the model is correct. Full 'WLS' is highly sensitive to model misspecification – it can heavily weight the fourth-order moments of the distribution, so small deviations between the observed fourth-order moments and those implied by the model can lead to poor estimates.

\strong{Behavior with all-continuous data}
When only continuous variables are present, the argument \code{allContinuousMethod} dictates how to process the data.

The default, \emph{cumulants} is a good choice for non-normal data. This uses the asymptotically distribution free (ADF) method of Browne (1984) and computes the fourth order \emph{cumulants} for the weight matrix: thus, the name. It is generally fast and ADF up to elliptical distributions. Data computed using cumulants should also be more accurate than via marginals (because the whole covariance is a single analytic expression, with no estimation involved). 

\emph{note}: The \emph{cumulants} method does not handle missing data. It also does not return weights or summary statistics for the means.

The alternative option, 'marginals', uses methods similar to those used in processing ordinal and joint ordinal-continuous data.  By contrast with cumulants, marginals returns weights and summary statistics for the means.

When data are not all continuous, \code{allContinuousMethod} is ignored, and means are modelled.

\emph{Usage Notes}:

Model results can be reported using the \link{summary} function, or accessed directly in the 'output' slot of the model (i.e., \code{model$output}). Components of the output may also be accessed and used in the same way, i.e., via the \code{$} and \code{[]} \code{\link{Extract}} functions.

Summary statistics are returned in the MxData object in an
\code{observedStats} list. If \code{observedStats} are already present
and in the appropriate shape then they are reused. It is also possible
to provide your own arbitrary user supplied observed statistics using
this same approach.
}


\value{
Returns a new MxFitFunctionWLS object. One and only one fit function object should be included in each model, along with an associated \link{mxExpectationNormal} or \code{\link{mxExpectationRAM}} object.
}

\references{
The OpenMx User's guide can be found at http://openmx.ssri.psu.edu/documentation.

Browne, M. W. (1984). Asymptotically distribution-free methods for the analysis of covariance structures. \emph{British Journal of Mathematical and Statistical Psychology}, \strong{37}, 62-83.
}

\seealso{
Other fit functions:
\code{\link{mxFitFunctionMultigroup}}, \code{\link{mxFitFunctionML}}, 
\code{\link{mxFitFunctionAlgebra}}, \code{\link{mxFitFunctionGREML}},
\code{\link{mxFitFunctionR}}, \code{\link{mxFitFunctionRow}}

More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}

\examples{

# Create and fit a WLS model using RAM, and then using matrices.

library(OpenMx)

# Simulate some data where y = .5x + error

x = rnorm(1000, mean = 0, sd = 1)
y = 0.5*x + rnorm(1000, mean = 0, sd = 1)
tmpFrame = data.frame(x, y)
varNames = names(tmpFrame)

# =======================
# = A RAM model example =
# =======================

m1 = mxModel("my_first_WLS", type="RAM",
	manifestVars = c("x", "y"),
	mxPath(c("x", "y"), arrows=2, values=1),
	mxPath("x", to="y"),
	mxFitFunctionWLS(),
	mxData(tmpFrame, 'raw')
)

m1 = mxRun(m1)
summary(m1)$parameters

originalSE <- m1$output$standardErrors

# Here are the summary statistics
print(m1$data$observedStats)

# Use a different weight matrix
os <- m1$data$observedStats
os$acov <- rWishart(1, nrow(tmpFrame)/2, diag(3))[,,1]
os$fullWeight <- os$acov
m1$data$observedStats <- os

# Can check whether our new weights are used
m1$data$verbose <- 1L

# Run model
m1 <- mxRun(m1)

# SE indeed changed due to new weights
print(m1$output$standardErrors - originalSE)

# ==========================
# = A matrix-based example =
# ==========================

# Define matrices for Symmetric (S) and Asymmetric (A) paths and an Identity matrix.

S <- mxMatrix(type = "Full", nrow = 2, ncol = 2, values=c(1,0,0,1), 
              free=c(TRUE,FALSE,FALSE,TRUE), labels=c("Vx", NA, NA, "Vy"), name = "S")
A <- mxMatrix(type = "Full", nrow = 2, ncol = 2, values=c(0,1,0,0), 
              free=c(FALSE,TRUE,FALSE,FALSE), labels=c(NA, "b", NA, NA), name = "A")
I <- mxMatrix(type="Iden", nrow=2, ncol=2, name="I")

# Build the model

tmpModel <- mxModel(model="exampleModel", 
	# Add the S, A, and I matrices constructed above
	S, A, I, 

	# Define the expectation
	mxAlgebra(name="expCov", solve(I-A) \%*\% S \%*\% t(solve(I-A))),

	# Choose a normal expectation and WLS as the fit function
	mxExpectationNormal(covariance= "expCov", dimnames= varNames),
	mxFitFunctionWLS(), 

   # Add the data
	mxData(tmpFrame, 'raw')
)

# Fit the model and print a summary
tmpModel <- mxRun(tmpModel)
summary(tmpModel)

}
