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

\name{mxBootstrapEval}
\alias{mxBootstrapEval}
\alias{omxBootstrapEvalCov}
\alias{omxBootstrapEval}
\alias{omxBootstrapEvalByName}
\alias{mxBootstrapEvalByName}

\title{Evaluate Values in a bootstrapped MxModel}

\description{
   This function can be used to evaluate an arbitrary R expression that includes named entities from a \link{MxModel} object, or labels from a \link{MxMatrix} object.
}

\usage{
mxBootstrapEval(expression, model, defvar.row = 1, ...,
 bq=c(.25,.75), method=c('bcbci','quantile'))

mxBootstrapEvalByName(name, model, defvar.row = 1, ...,
 bq=c(.25,.75), method=c('bcbci','quantile'))

omxBootstrapEval(expression, model, defvar.row = 1L, ...)

omxBootstrapEvalCov(expression, model, defvar.row = 1L, ...)

omxBootstrapEvalByName(name, model, defvar.row=1L, ...)
}

\arguments{
   \item{expression}{An arbitrary R expression.}
   \item{name}{The character name of an object to evaluate.}
   \item{model}{The model in which to evaluate the expression.}
   \item{defvar.row}{The row to use for definition variables when compute=TRUE (defaults to 1).  When compute=FALSE, values for definition variables are always taken from the first (i.e., first before any automated sorting is done) row of the raw data.}
  \item{...}{Not used.  Forces remaining arguments to be specified by name.}
\item{bq}{numeric. A vector of bootstrap quantiles at which to summarize the bootstrap replication.}
\item{method}{character. One of \sQuote{quantile} or \sQuote{bcbci}.}
}

\details{
The argument \sQuote{expression} is an arbitrary R expression.  Any named entities that are used within the R expression are translated into their current value from the model. Any labels from the matrices within the model are translated into their current value from the model. Finally the expression is evaluated and the result is returned.  To enable debugging, the \sQuote{show} argument has been provided.  The most common mistake when using this function is to include named entities in the model that are identical to R function names.  For example, if a model contains a named entity named \sQuote{c}, then the following mxEval call will return an error: \code{mxEval(c(A, B, C), model)}.

The \code{mxEvalByName} function is a wrapper around \code{mxEval} that takes a character instead of an R expression.

\emph{nb}: \sQuote{bcbci} stands for \sQuote{bias-corrected bootstrap confidence interval}

The default behavior is to use the \sQuote{bcbci} \code{method}, due to its superior theoretical properties.

}

\value{

\code{omxBootstrapEval} and \code{omxBootstrapEvalByName} return the raw matrix of
\code{cvectorize}'d results. \code{omxBootstrapEvalCov} returns the
covariance matrix of the \code{cvectorize}'d results.
\code{mxBootstrapEval} and \code{mxBootstrapEvalByName} return
the \code{cvectorize}'d results summarized by \code{method} at quantiles \code{bq}.
}

\references{
The OpenMx User's guide can be found at \url{https://openmx.ssri.psu.edu/documentation/}.
}

\seealso{
\link{mxAlgebra} to create algebraic expressions inside your model 
and \link{mxModel} for the model object mxEval looks inside when
evaluating. \link{mxBootstrap} to create bootstrap data.
}

\examples{

library(OpenMx)
# make a unit-weighted 10-row data set of values 1 thru 10
myData = mxData(data.frame(weight=1.0, value=1:10), "raw", weight = "weight")
sum(1:10)

# Model sums data$value (sum(1:10)= 55), subtracts "A", squares the result, 
# and tries to minimize this (achieved by setting A=55)
testModel = mxModel(model = "testModel1", myData,
	mxMatrix(name  = "A", "Full", nrow = 1, ncol = 1, values = 1, free=TRUE),
	# nb: filteredDataRow is an auto-generated matrix of
	# non-missing data from the present row.
	# This is placed into the "rowResults" matrix (also auto-generated)
	mxAlgebra(name = "rowAlg", data.weight * filteredDataRow),
	# Algebra to turn the rowResults into a single number
	mxAlgebra(name = "reduceAlg", (sum(rowResults) - A)^2),
	mxFitFunctionRow(
		rowAlgebra    = "rowAlg",
		reduceAlgebra = "reduceAlg",
		dimnames      = "value"
	)
	# no need for an MxExpectation object when using mxFitFunctionRow
)

testModel = mxRun(testModel) # A is estimated at 55, with SE= 1
testBoot = mxBootstrap(testModel)
summary(testBoot) # A is estimated at 55, with SE= 0

# Let's compute A^2 (55^2 = 3025)
mxBootstrapEval(A^2, testBoot)
#      SE 25.0% 75.0%
# [1,]  0  3025  3025

}
