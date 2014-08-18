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

\name{mxRestore}
\alias{mxRestore}

\title{Restore From Checkpoint File}

\description{
   The function loads the last saved state from a checkpoint file.
}

\usage{
mxRestore(model, chkpt.directory = ".", chkpt.prefix = "")
}

\arguments{
    \item{model}{\link{MxModel} object to be loaded.}    
    \item{chkpt.directory}{character. Directory where the checkpoint file is located.}
    \item{chkpt.prefix}{character. Prefix of the checkpoint file.}
}

\details{
In general, the arguments \sQuote{chkpt.directory} and \sQuote{chkpt.prefix} should be identical to the \code{\link{mxOption}}: \sQuote{Checkpoint Directory} and \sQuote{Checkpoint Prefix} that were specificed on the model before execution.

Alternatively, the checkpoint file can be manually loaded as a data.frame in R.  Use \code{\link{read.table}} with the options \sQuote{header=TRUE}, \sQuote{sep="\t"}, \sQuote{stringsAsFactors=FALSE} and \sQuote{check.names=FALSE}.
}

\value{
   Returns an MxModel object with free parameters updated to the last saved values.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{

library(OpenMx)

# Simulate some data

x=rnorm(1000, mean=0, sd=1)
y= 0.5*x + rnorm(1000, mean=0, sd=1)
tmpFrame <- data.frame(x, y)
tmpNames <- names(tmpFrame)

# Create a model that includes an expected covariance matrix,
# an expectation function, a fit function, and an observed covariance matrix

data <- mxData(cov(tmpFrame), type="cov", numObs = 1000)
expCov <- mxMatrix(type="Symm", nrow=2, ncol=2, values=c(.2,.1,.2), free=TRUE, name="expCov")
expFunction <- mxExpectationNormal(covariance="expCov", dimnames=tmpNames)
fitFunction <- mxFitFunctionML()
testModel <- mxModel(model="testModel", expCov, data, expFunction, fitFunction)

#Use mxRun to optimize the free parameters in the expected covariance matrix
modelOut <- mxRun(testModel, checkpoint = TRUE)
modelOut$expCov

#Use mxRestore to load the last checkpoint saved state of the model
modelRestore <- mxRestore(testModel)
modelRestore$expCov


}
