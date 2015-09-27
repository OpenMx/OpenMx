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

\name{mxGetExpected}
\alias{mxGetExpected}
\alias{imxGetExpectationComponent}


\title{Extract the component from a model's expectation}

\description{
This function extracts the expected means, covariance, or thresholds from a model.
}

\usage{
mxGetExpected(model, component, defvar.row=1)
imxGetExpectationComponent(model, component, defvar.row=1)
}


\arguments{
   \item{model}{MxModel object from which to extract the expectation component.}
   \item{component}{A character. The name of the component to extract.}
   \item{defvar.row}{A row index. Which row to load for definition variables.}
}

\details{
The expected means, covariance, or thresholds can be extracted from
Normal (\link{mxExpectationNormal}), RAM (\link{mxExpectationRAM}), and
LISREL (\link{mxExpectationLISREL}) models. When more than one component
is requested, the components will be returned as a list.

If component 'vector' is requested then the non-redundent coefficients
of the expected manifest distribution will be returned as a vector.
}

\value{
See details.
}


\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{
   
#----------
require(OpenMx)
manifests <- paste("x", 1:5, sep="")
latents <- c("G")
factorModel <- mxModel("One Factor",
      type="RAM",
      manifestVars = manifests,
      latentVars = latents,
      mxPath(from=latents, to=manifests),
      mxPath(from=manifests, arrows=2),
      mxPath(from=latents, arrows=2,
            free=FALSE, values=1.0),
      mxPath(from = 'one', to = manifests))

mxGetExpected(factorModel, "covariance")
# oops.  Starting values indicate a zero covariance matrix.
#  Probably should adjust them.

}
