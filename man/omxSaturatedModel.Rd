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

\name{omxSaturatedModel}
\alias{omxSaturatedModel}
\alias{mxRefModels}

\title{Create Reference (Saturated and Independence) Models}

\description{
   This function creates and optionally runs saturated and independence (null) models of a base model or data set for use with \link{mxSummary} to obtain more fit indices.
}


\usage{
mxRefModels(x, run=FALSE)
}

\arguments{
   \item{x}{A MxModel object, data frame, or matrix.}
   \item{run}{logical. If TRUE runs the models before returning; otherwise returns built models without running.}
}

\details{
For typical structural equation models the saturated model is the free-est possible model.  All covariances and, when possilbe, all means are estimated.  In the case of ordinal data, the ordinal means are fixed to zero and the thresholds are estimated.  When the ordinal data are binary, those variances are also constrained to one.  This is the free-est possible model, only constrained for model identification.  The saturated model is used to create the RMSEA, and Chi-squared fit indices.

The independence model, sometimes called the null model, is a model of each variable being completely independent of every other variable.  As such, all the variances and, when possible, all means are estimated.  However, covariances are set to zero.  Ordinal variables are handled the same for the independence and saturated models.  The independence model is used, along with the saturated model, to create CFI and TLI fit indices.

When the \link{mxFitFunctionMultigroup} fit function is used, \code{mxRefModels} creates the appropriate multigroup saturated and independence models.  Saturated and independence models are created separately for each group.  Each group has its own saturated and independence model.  The multigroup saturated model is a multigroup model where each group has its own saturated model, and similarly for the independence model.

One potentially important limitation of the \code{mxRefModels} function is for behavioral genetics models.  If variables 'x', 'y', and 'z' are measured on twins 1 and 2 creating the modeled variables 'x1', 'y1', 'z1', 'x2', 'y2', 'z2', then this function may not create the intended saturated or independence models.  In particular, the means of 'x1' and 'x2' are estimated separately.  Similarly, the covariance of 'x1' with 'y1' and 'x2' with 'y2' are allowed be be distinct: \eqn{cov(x1, y1) != cov{x2, y2}}.  Moreover, the cross-twin covariances are estimated: e.g. \eqn{cov(x1, y2) != 0}.

Another potential misuse of this function is for models with definition variables.  If definition variables are used, the saturated and independence model may not be correct because they do not account for the definition variables.

When an MxModel has been run, some effort is made to make the reference models for only the variables used in the model.  For covariance data, all variables are modeled by default.  For raw data when the model has been run, only the modeled variables are used in the reference models.  This matches the behavior of \link{mxModel}.

In general, it is best practice to give \code{mxRefModels} a model that has already been run.

For IFA models (\link{mxExpectationBA81}), the independence model
preserves equality constraints among item parameters from the original model.

\code{mxRefModels()} is not compatible with \link[OpenMx:mxExpectationGREML]{GREML expectation}, as there is no sensible general definition for a saturated GREML-type model.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
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
                       mxData(observed=cov(demoOneFactor), type="cov", numObs=500))
summary(factorRun <- mxRun(factorModel))
factorSat <- mxRefModels(factorRun, run=TRUE)
summary(factorRun, refModels=factorSat)
}

