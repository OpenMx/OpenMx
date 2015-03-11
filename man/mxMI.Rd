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

\name{mxMI}
\alias{mxMI}


\title{Estimate Modification Indices for MxModel Objects}

\description{
This function estimates the change in fit function value resulting from freeing currently fixed parameters.
}

\usage{
mxMI(model, matrices=NA, full=TRUE)
}


\arguments{
   \item{model}{An MxModel for which modification indices are desired.}
   \item{matrices}{Character vector.  The names of the matrices in which to search for modification}
   \item{full}{Logical.  Whether or not to return the full modification index in addition to the restricted.}
}

\details{
Modification indices provide an estimate of how much the fit function value would change if a parameter that is currently fixed was instead freely estimated.  There are two versions of this estimate: a restricted version and an full version.  The restricted version is reported as the MI and is much faster to compute.  The full version is reported as MI.Full.  The full version accounts for the \emph{total} change in fit function value resulting from the newly freed parameter.  The restricted version only accounts for the change in the fit function due to the movement of the new free parameter.  In particular, the restricted version does not account for the change in fit function value due to the other free parameters moving in response to the new parameter.

When the fit function is in minus two log likelihood units (e.g. \link{\code{mxFitFunctionML}}), then the MI will be approximately chi squared distributed with 1 degree of freedom.  Using a p-value of 0.01 has been suggested.  Hence, a MI greater than \code{qchisq(p=1-0.01, df=1)}, or 6.63, is suggestive of a modification.

Users should be cautious in their use of modification indices.  If a model was created with the aid of MIs, then it should \emph{always} be reported.  \emph{Do not pretend that you have a theoretical reason for part of a model that was put there because it was suggested by a modification index.  This is fraud.}  When using modification indices there are two options for best practices.  First, you can report the analyses as exploratory.  Document all the explorations that you did, and know that your results may or may not generalize.  Second, you can use cross-validation.  Reserve part of your data for exploration, and use the remaining data to test if the exploratory model generalizes to new data.
}

\value{
    A named list with components
    \describe{
        \item{MI}{The restricted modification index.}
        \item{MI.Full}{The full modification index.}
        \item{plusOneParamModels}{A list of models with one additional free parameter}
    }
}


\references{
S\ifelse{latex}{\out{\"o}}{\ifelse{html}{\out{&ouml;}}{o}}rbom, D.  (1989).  Model Modification.  \emph{Psychometrika, 54}, 371-384.

The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{
# Create a model
require(OpenMx)
data(demoOneFactor)
manifests <- names(demoOneFactor)
latents <- c("G")
factorModel <- mxModel("One Factor",
      type="RAM",
      manifestVars = manifests,
      latentVars = latents,
      mxPath(from=latents, to=manifests),
      mxPath(from=manifests, arrows=2),
      mxPath(from=latents, arrows=2,
            free=FALSE, values=1.0),
      mxPath(from = 'one', to = manifests),
      mxData(demoOneFactor, type="raw"))
factorRun <- mxRun(factorModel)

# See if it should be modified
fim <- mxMI(factorRun, full=FALSE)
round(fim$MI, 3)
plot(fim$MI, ylim=c(0, 10))
abline(h=qchisq(p=1-0.01, df=1)) # line of "significance"

}
