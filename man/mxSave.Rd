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

\name{mxSave}
\alias{mxSave}

\title{Save End State to Checkpoint File}

\description{
   The function saves the last state of a model to a checkpoint file.
}

\usage{
mxSave(model, chkpt.directory = ".", chkpt.prefix = "")
}

\arguments{
    \item{model}{\link{MxModel} object to be loaded.}    
    \item{chkpt.directory}{character. Directory where the checkpoint file is located.}
    \item{chkpt.prefix}{character. Prefix of the checkpoint file.}
}

\details{
In general, the arguments \sQuote{chkpt.directory} and \sQuote{chkpt.prefix} should be identical to the \code{\link{mxOption}}: \sQuote{Checkpoint Directory} and \sQuote{Checkpoint Prefix} that were specificed on the model before execution.

Alternatively, the checkpoint file can be manually loaded as a data.frame in R.  Use \code{\link{read.table}} with the options \code{header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE}.
}

\value{
   Returns a logical indicating the succes of writing the checkpoint file to the checkpoint directory.
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
modelOut <- mxRun(testModel)
modelOut$expCov

# Save the ending state of modelOut in a checkpoint file
mxSave(modelOut)

# Restore the saved model from the checkpoint file
modelSaved <- mxRestore(testModel)
modelSaved$expCov

# Imagine how much time you saved by not having to re-run the
# model that took hours or days to run.
}
