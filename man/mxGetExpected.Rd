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

\name{mxGetExpected}
\alias{mxGetExpected}
\alias{imxGetExpectationComponent}


\title{Extract the component from a model's expectation}

\description{
This function extracts the expected means, covariance, or thresholds from a model.
}

\usage{
mxGetExpected(model, component, defvar.row=1, subname=model$name)
imxGetExpectationComponent(model, component, defvar.row=1, subname=model$name)
}


\arguments{
   \item{model}{MxModel object from which to extract the expectation component.}
   \item{component}{Character vector. The name(s) of the component(s) to extract.  Recognized names are \dQuote{covariance}, \dQuote{means}, and \dQuote{thresholds}.}
   \item{defvar.row}{A row index. Which row to load for definition variables.}
   \item{subname}{Name of the submodel to evaluate.}
}

\details{
The expected means, covariance, or thresholds can be extracted from
Normal (\link{mxExpectationNormal}), RAM (\link{mxExpectationRAM}), and
LISREL (\link{mxExpectationLISREL}) models. When more than one component
is requested, the components will be returned as a list.

If component 'vector' is requested then the non-redundant coefficients
of the expected manifest distribution will be returned as a vector.

If component 'standVector' is requested then the same parameter structure as
'vector' is returned, but it is standardized.  For Normal expectations the covariances
are returned as correlations, the means are returned as zeros, and the
thresholds are returned as z-scores.  For the thresholds the z-scores
are computed by using the model-implied means and variances.

Note that capitalization is ignored for the 'standVector' option, so 'standvector'
is also acceptable.
}

\value{
See details.
}


\references{
The OpenMx User's guide can be found at http://openmx.ssri.psu.edu/documentation.
}

\examples{
   
# ===============================================
# = Build a 1-factor CFA, with bad start values =
# ===============================================
require(OpenMx)
manifests = paste("x", 1:5, sep="")
latents = c("G")
factorModel = mxModel("One Factor", type="RAM",
      manifestVars = manifests,
      latentVars = latents,
      mxPath(from = latents, to = manifests),
      mxPath(from = manifests, arrows = 2),
      mxPath(from = latents, arrows = 2, free = FALSE, values = 1.0),
      mxPath(from = 'one', to = manifests),
	  mxData(demoOneFactor, type = "raw")
)

# ============================================================================
# = What do our starting values indicate about the expected data covariance? =
# ============================================================================
mxGetExpected(factorModel, "covariance")

# Oops. Starting values indicate an expected zero-covariance matrix.
# The model likely won't run from these start values.
# Let's adjust them:

factorModel = mxModel("One Factor", type = "RAM",
      manifestVars = manifests, latentVars = latents,
      # Reasonable start VALUES
	  mxPath(from = latents, to = manifests, values = .2),
      mxPath(from = manifests, arrows = 2),
      mxPath(from = latents, arrows = 2, free = FALSE, values = 1.0),
      mxPath(from = 'one', to = manifests),
	  mxData(demoOneFactor, type = "raw")
)

mxGetExpected(factorModel, "covariance")
#      x1   x2   x3   x4   x5
# x1 0.04 0.04 0.04 0.04 0.04
# x2 0.04 0.04 0.04 0.04 0.04
# x3 0.04 0.04 0.04 0.04 0.04
# x4 0.04 0.04 0.04 0.04 0.04
# x5 0.04 0.04 0.04 0.04 0.04

# And this version will run:
factorModel = mxRun(factorModel)

}
