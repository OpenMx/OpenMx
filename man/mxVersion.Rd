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

\name{mxVersion}
\alias{mxVersion}

\title{Returns Current Version String}

\description{
   This function returns a string with the current version number of OpenMx.
   Optionally (with verbose = TRUE (the default)), it prints a message containing
   the version of R, the platform, and the optimiser.
}

\usage{
mxVersion(model = NULL, verbose = TRUE)
}

\arguments{
  \item{model}{ optional \code{\link{MxModel}} to request optimizer from (default = NULL)}
  \item{verbose}{Whether to print version information to the console (default = TRUE)}
}
\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{

# Print useful version information.
mxVersion()
# If you just want the version, use this call.
x = mxVersion(verbose=FALSE)

library(OpenMx)
data(demoOneFactor)  # load the demoOneFactor dataframe
manifests <- names(demoOneFactor) # set the manifest to the 5 demo variables
latents   <- c("G")  # define 1 latent variable
model <- mxModel(model = "One Factor", type = "RAM",
    manifestVars = manifests,
    latentVars   = latents,
    mxPath(from = latents, to = manifests, labels = paste("b", 1:5, sep = "")), 
    mxPath(from = manifests, arrows = 2  , labels = paste("u", 1:5, sep = "")), 
    mxPath(from = latents  , arrows = 2  , free = FALSE, values = 1.0), 
    mxData(cov(demoOneFactor), type = "cov", numObs = 500)
)
mxVersion(model, verbose = TRUE)

}
