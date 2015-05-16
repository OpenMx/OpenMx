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

\name{mxRun}
\alias{mxRun}

\title{Send a Model to the Optimizer}

\description{
   This function begins optimization on the top-level model.
}

\usage{
mxRun(model, ..., intervals = NULL, silent = FALSE, suppressWarnings = FALSE, 
    unsafe = FALSE, checkpoint = FALSE, useSocket = FALSE, onlyFrontend = FALSE, 
    useOptimizer = TRUE)
}

\arguments{
    \item{model}{A \link{MxModel} object to be optimized.}    
    \item{...}{Not used.  Forces remaining arguments to be specified by name.}
    \item{intervals}{A boolean indicating whether to compute the specified confidence intervals.}
    \item{silent}{A boolean indicating whether to print status to terminal.}
    \item{suppressWarnings}{A boolean indicating whether to suppress warnings.}    
    \item{unsafe}{A boolean indicating whether to ignore errors.}
    \item{checkpoint}{A boolean indicating whether to periodically write parameter values to a file.}
    \item{useSocket}{A boolean indicating whether to periodically write parameter values to a socket.}
    \item{onlyFrontend}{A boolean indicating whether to run only front-end model transformations.}
    \item{useOptimizer}{A boolean indicating whether to run only the log-likelihood of the current free parameter values but not move any of the free parameters.}
}

\details{
The mxRun function is used to optimize free parameters in \link{MxModel} objects based on an expectation function and fit function. MxModel objects included in the mxRun function must include an appropriate expectation and fit functions.

If the \sQuote{silent} flag is TRUE, then model execution will not print any status messages to the terminal.

If the \sQuote{suppressWarnings} flag is TRUE, then model execution will not issue a warning if NPSOL returns a non-zero status code.

If the \sQuote{unsafe} flag is TRUE, then any error conditions will throw a warning instead of an error.  It is strongly recommended to use this feature only for debugging purposes.

Free parameters are estimated or updated based on the expectation and fit functions. These estimated values, along with estimation information and model fit, can be found in the 'output' slot of MxModel objects after mxRun has been used.

If a model is dependent on or shares parameters with another model, both models must be included as arguments in another MxModel object. This top-level MxModel object must include expectation and fit functions in both submodels, as well as an additional fit function describing how the results of the first two should be combined.
}

\value{
   Returns an MxModel object with free parameters updated to their final values.
   The return value contains an "output" slot with the results of 
   optimization.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{
# Create and run the 1-factor CFA on the openmx.psyc.virginia.edu front page

library(OpenMx)

data(demoOneFactor)  # load the demoOneFactor dataframe

manifests <- names(demoOneFactor) # set the manifest to the 5 demo variables
latents   <- c("G")  # define 1 latent variable
model <- mxModel(model="One Factor", type="RAM",
    manifestVars = manifests,
    latentVars   = latents,
    mxPath(from=latents, to=manifests, labels=paste("b", 1:5, sep="")),
    mxPath(from=manifests, arrows=2, labels=paste("u", 1:5, sep="")),
    mxPath(from=latents  , arrows=2, free=FALSE, values=1.0),
    mxData(cov(demoOneFactor), type="cov", numObs=500)
)
model <- mxRun(model) # Run the model, returning the result into model
summary(model) # Show summary of the fitted model

}
