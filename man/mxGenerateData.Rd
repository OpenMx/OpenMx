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

\name{mxGenerateData}
\alias{mxGenerateData}


\title{Generate data based on an mxModel (or a data.frame)}

\description{
This function returns a new (simulated) data set based on either the model-implied distribution if
a model is provided, OR saturated model if a data.frame is given in the model parameter.

See below for important details
}

\usage{
mxGenerateData(model, nrows, returnModel=FALSE, use.miss = TRUE,
 ..., .backend=TRUE, subname=NULL, empirical=FALSE, nrowsProportion)
}


\arguments{
   \item{model}{A data.frame or MxModel object upon which the data are generated.}
   \item{nrows}{Numeric. The number of rows of data to generate (default = same as in the original data)}
   \item{returnModel}{ Whether to return the model with new data, or just return the new data.frames (default)}
   \item{use.miss}{ Whether to approximate the missingness pattern of the original data (TRUE by default).}
\item{...}{Not used; forces remaining arguments to be specified by name.}
\item{.backend}{Whether to use the backend to generate data (TRUE by default for speed)}
\item{subname}{If given, limits data generation to this sub model.}
\item{empirical}{Whether the generate data should match the distribution
   of the current data exactly. Uses \link[MASS]{mvrnorm} instead of
   \link[mvtnorm]{rmvnorm}}
 \item{nrowsProportion}{Numeric. The number of rows of data to generate
   expressed as a proportion of the current number of rows.}
}

\details{
When given a data.frame as a model, the model is assumed to be
saturated multivariate Gaussian and the expected distribution
is obtained using \link{mxDataWLS}.  In this case, the default number of rows is assumed to be the number of rows in the original data.frame, but any other number of rows can also be requested.

When given an MxModel, 
the model-implied means and covariance are extracted.  It then generates data with the same mean and covariance.  Data can be generated based on Normal (\link{mxExpectationNormal}), RAM (\link{mxExpectationRAM}), LISREL (\link{mxExpectationLISREL}), and state space (\link{mxExpectationStateSpace}) models.

Please note that this function samples data from the model-implied distribution(s); it does not sample from the data object in the model.  That is, this function generates new data rather than pulling data that already exist from the model.

Thresholds and ordinal data are implemented by generating continuous data and then using \link{cut} and \link{mxFactor} to break the continuous data at the thresholds into an ordered factor.

If the model has definition variables, then a data set must be included in the model object and the number of rows requested must match the number of rows in the model data.  In this case the means, covariance, and thresholds are reevaluated for each row of data, potentially creating a a different mean, covariance, and threshold structure for every generated row of data.

For state space models (i.e. models with an \link{mxExpectationStateSpace} or \link{mxExpectationStateSpaceContinuousTime} expectation), the data are generated based on the autoregressive structure of the model.  The rows of data in a state space model are not independent replicates of a stationary process.  Rather, they are the result of a latent (possibly non-stationary) autoregressive process.  For state space models different rows of data often correspond to different times.  As alluded to above, data generation works for discrete time state space models and hybrid continuous-discrete time state space models.  The latter have a continuous process that is measured as discrete times.

The \code{subname} parameter is used to limit data generation to
the given submodel. The reason you wouldn't pass the submodel in
the \code{model} argument is that some parts of the submodel might
depend on objects in other submodels that are part of the model.
}

\value{
    A data.frame, list of data.frames, or model populated with the new data
    (depending on the \code{returnModel} parameter).
    Raw data is always returned even if the original model contained
    covariance or some other non-raw data.
}


\references{
The OpenMx User's guide can be found at http://openmx.ssri.psu.edu/documentation.
}

\examples{

# ====================================
# = Demonstration for empirical=TRUE =
# ====================================
popCov <- cov(Bollen[, 1:8])*(nrow(Bollen)-1)/nrow(Bollen)
got <- mxGenerateData(Bollen[, 1:8], nrows=nrow(Bollen), empirical = TRUE)
cov(got) - popCov  # pretty close, given 8 variables to juggle!
round(cov2cor(cov(got)) - cov2cor(popCov), 4)

# ===========================================
# = Create data based on state space model. =
# ===========================================

require(OpenMx)
nvar <- 5
varnames <- paste("x", 1:nvar, sep="")
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

# Add simulated data to model and run
ssModel <- mxModel(ssModel, mxData(ssData, 'raw'))
ssRun <- mxRun(ssModel)

# Compare parameters from random data to the generating model
cbind(Rand = omxGetParameters(ssRun), Gen = omxGetParameters(ssModel))

# Note the parameters should be "close" (up to sampling error)
# to the generating values


# =========================================
# = Demo generating new data from a model =
# =========================================
require(OpenMx)
manifests <- paste0("x", 1:5)
originalModel <- mxModel("One Factor", type="RAM",
      manifestVars = manifests,
      latentVars = "G",
      mxPath(from="G", to=manifests, values=.8),
      mxPath(from=manifests, arrows=2, values=.2),
      mxPath(from="G"  , arrows=2, free=FALSE, values=1.0),
      mxPath(from = 'one', to = manifests)
)

factorData <- mxGenerateData(originalModel, 1000)
newData = mxData(cov(factorData), type="cov", 
	numObs=nrow(factorData), means = colMeans(factorData)
)
newModel <- mxModel(originalModel, newData)
newModel <- mxRun(newModel)
cbind(
	Original = omxGetParameters(originalModel),
	Generated = round(omxGetParameters(newModel), 4),
	Delta = round(
		omxGetParameters(originalModel) - 
		omxGetParameters(newModel), 3)
)

# And again with empirical = TRUE

factorData <- mxGenerateData(originalModel, 1000, empirical = TRUE)
newData = mxData(cov(factorData),
	type = "cov", 
	numObs = nrow(factorData),
	means = colMeans(factorData)
)

newModel <- mxModel(originalModel, newData)
newModel <- mxRun(newModel)

cbind(
	Original  = omxGetParameters(originalModel),
	Generated = round(omxGetParameters(newModel), 4),
	Delta     = omxGetParameters(originalModel) - 
		 	    omxGetParameters(newModel)
)

}
