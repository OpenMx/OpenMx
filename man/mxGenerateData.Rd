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

\name{mxGenerateData}
\alias{mxGenerateData}


\title{Generate data based on an MxModel object}

\description{
This function creates a randomly sampled data set based on the model.
}

\usage{
mxGenerateData(model, nrows)
}


\arguments{
   \item{model}{An MxModel object upon which the data are generated.}
   \item{nrows}{Numeric. The number of rows of data to generate.}
}

\details{
This function looks inside the MxModel object to extract the model-implied means and covariance.  It then generates data with the same mean and covariance.  Data can be generated based on Normal (\link{mxExpectationNormal}), RAM (\link{mxExpectationRAM}), and LISREL (\link{mxExpectationLISREL}) models.

Thresholds and ordinal data are not yet implemented.

For state space models (i.e. models with an \link{mxExpectationStateSpace} expectation), the data are generated based on the autoregressive structure of the model.  The rows of data in a state space model are not independent replicates of a stationary process.  Rather, they are the result of a latent (possibly non-stationary) autoregressive process.  For state space models different rows of data often correspond to different times.
}

\value{
    A matrix or data.frame with nrows rows.
}


\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{

#----------
# Create data based on state space model.
require(OpenMx)
varnames <- paste("x", 1:5, sep="")
ssModel <- mxModel(model="State Space Manual Example",
    mxMatrix("Full", 1, 1, TRUE, .3, name="A"),
    mxMatrix("Zero", 1, 1, name="B"),
    mxMatrix("Full", nvar, 1, TRUE, .6, name="C", dimnames=list(varnames, "F1")),
    mxMatrix("Zero", nvar, 1, name="D"),
    mxMatrix("Diag", 1, 1, FALSE, 1, name="Q"),
    mxMatrix("Diag", nvar, nvar, TRUE, .2, name="R"),
    mxMatrix("Zero", 1, 1, name="x0"),
    mxMatrix("Diag", 1, 1, FALSE, 1, name="P0"),
    mxMatrix("Zero", 1, 1, name="u"),
    mxExpectationStateSpace("A", "B", "C", "D", "Q", "R", "x0", "P0", "u"),
    mxFitFunctionML()
)

ssData <- mxGenerateData(ssModel, 200) # 200 time points

# Add simulated data to model
ssModel <- mxModel(ssModel, mxData(ssData, 'raw'))

# Fit model to simulated data
ssRun <- mxRun(ssModel)

# Compare parameters estimated from random data to
#  their true generating values
cbind(Rand=omxGetParameters(ssRun), Gen=omxGetParameters(ssModel))
# Note the parameters should be "close" (up to sampling error)
# to the generating values


#----------
require(OpenMx)
manifests <- paste("x", 1:5, sep="")
latents <- c("G")
factorModel <- mxModel("One Factor",
      type="RAM",
      manifestVars = manifests,
      latentVars = latents,
      mxPath(from=latents, to=manifests, values=.8),
      mxPath(from=manifests, arrows=2, values=.2),
      mxPath(from=latents, arrows=2,
            free=FALSE, values=1.0),
      mxPath(from = 'one', to = manifests))

factorData <- mxGenerateData(factorModel, 100)

factorModel <- mxModel(factorModel, mxData(factorData, 'raw'))
factorRun <- mxRun(factorModel)
cbind(Rand=omxGetParameters(factorRun), Gen=omxGetParameters(factorModel))


}
