%
%   Copyright 2007-2017 The OpenMx Project
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

\name{omxParallelCI}
\alias{omxParallelCI}
\alias{omxRunCI}

\title{
Dedicated functions for calculating confidence intervals without redoing the primary optimization.
}

\description{
\code{omxParallelCI()} calculates confidence intervals for an already-\link[=mxRun]{run} \link{MxModel} object that contains an \link{MxInterval} object (i.e., an \code{\link{mxCI}()} statement), without recalculating point estimates, fitfunction derivatives, or expectations.  The "parallel" in the functions name refers to the not-yet-implemented feature of running independent \link[=mxModel]{submodels} in parallel.  \code{omxRunCI()} is a wrapper for \code{omxParallelCI()} with arguments \code{run=TRUE} and \code{independentSubmodels=FALSE}, and is the recommended interface.
}
\usage{
omxParallelCI(model, run = TRUE, verbose=0L, independentSubmodels=TRUE)

omxRunCI(model, verbose=0L)
}

\arguments{
  \item{model}{An \link{MxModel} object that contains an \link{MxInterval} object (i.e., an \code{\link{mxCI}()} statement).}
  \item{run}{Logical; should the resulting MxModel object be passed through \code{\link{mxRun}()} before being returned?  Defaults to \code{TRUE} (for \code{omxParallelCI()}).}
  \item{verbose}{Integer; defaults to zero; verbosity level passed to MxCompute* objects.}
  \item{independentSubmodels}{Logical; defaults to \code{TRUE} (for \code{omxParallelCI()}).  See below under "Details."}
}
\details{
When \code{independentSubmodels=TRUE}, \code{omxParallelCI()} creates an independent MxModel object for each quantity specified in the 'reference' slot of \code{model}'s \link{MxInterval} object, and places these independent MxModels inside \code{model}.  Each of these independent submodels calculates the confidence limits of its own quantity when the container is run.  When \code{independentSubmodels=FALSE}, no submodels are added to \code{model}.  Instead, \code{model} is provided with a dedicated \link[=mxComputeSequence]{compute plan} consisting only of an \link[=mxComputeConfidenceInterval][MxComputeConfidenceInterval] step.  Note that using \code{independentSubmodels=FALSE} will overwrite any compute plan already inside \code{model}.
}
\value{
The functions return \code{model}, augmented with independent submodels (if \code{independentSubmodels=TRUE}) or with a non-default compute plan (if \code{independentSubmodels=FALSE}), and possibly having been passed through \code{\link{mxRun}()} (if \code{run=TRUE}).  Naturally, if \code{run=FALSE}, the user can subsequently \link[=mxRun]{run} the returned model to obtain confidence intervals.  Users are cautioned that the returned model may not be very amenable to being further modified and re-fitted (e.g., having some free parameters fixed via \code{\link{omxSetParameters}()} and passed through \code{\link{mxRun}()} to get new point estimates) unless the added submodels or the non-default compute plan are eliminated.  The exception is if \code{run=TRUE} and \code{independentSubmodels=TRUE} (which is always the case with \code{omxRunCI()}), since the non-default compute plan is set to be non-persistent, and will automatically be replaced with a default compute plan the next time the model is passed to \code{\link{mxRun}()}.
}

\seealso{
\code{\link{mxCI}()}, \link{MxInterval}, \code{\link{mxComputeConfidenceInterval}}
}
\examples{
require(OpenMx)
data(demoOneFactor)
manifests <- names(demoOneFactor)
latents <- c("G")
factorModel <- mxModel("One Factor", 
	type="RAM",
	manifestVars=manifests, 
	latentVars=latents,
	mxPath(from=latents, to=manifests),
 	mxPath(from=manifests, arrows=2),
	mxPath(from=latents, arrows=2, free=FALSE, values=1.0),
	mxData(observed=cov(demoOneFactor), type="cov",
		numObs=500),
# add confidence intervals for free params in A and S matrices
	mxCI(c('A', 'S')))
factorRun <- mxRun(factorModel)
factorCI <- omxParallelCI(factorRun) 
factorCI2 <- omxRunCI(factorRun)
}
