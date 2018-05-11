%
%   Copyright 2007-2018 by the individuals mentioned in the source code history
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

\name{omxModelDeleteData}
\alias{omxModelDeleteData}

\title{Remove all instances of data from a model}

\description{
For very large data, it can be desirable to discard data after the model
is run. That is what the purpose of this function.

Data is discarded from the model and all submodels recursively.
}

\usage{
omxModelDeleteData(model)
}

\arguments{
   \item{model}{a MxModel object.}
}

\examples{

library(OpenMx)

data(demoOneFactor)
manifests <- names(demoOneFactor)
latents   <- c("G")
factorModel <- mxModel(model="One Factor", type="RAM",
      manifestVars = manifests,
      latentVars   = latents,
      mxPath(from=latents, to=manifests),
      mxPath(from=manifests, arrows=2),
      mxPath(from=latents, arrows=2,free=FALSE, values=1.0),
      mxData(cov(demoOneFactor), type="cov",numObs=500)
)
factorFit <-mxRun(factorModel)
object.size(factorFit)
factorFit <- omxModelDeleteData(factorFit)
object.size(factorFit)
factorFit$data

}
