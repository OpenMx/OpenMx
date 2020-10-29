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

\name{omxSaturatedModel}
\alias{omxSaturatedModel}
\alias{mxRefModels}

\title{Create Reference (Saturated and Independence) Models}

\description{
   This function creates and, optionally, runs saturated and independence (null) models of a base model or data set for use with \link{mxSummary} to enable fit indices that depend on these models.  Note that there are cases where this function is not valid for use, or should be used with caution (see below, under "Warnings").
}


\usage{
mxRefModels(x, run=FALSE, ..., distribution="default", equateThresholds = TRUE)
}

\arguments{
   \item{x}{A MxModel object, data frame, or matrix.}
   \item{run}{logical. If TRUE, runs the models before returning;
     otherwise returns built models without running.}
   \item{...}{Not used.  Forces remaining arguments to be specified by name.}
   \item{distribution}{character. Which distribution to assume.}
   \item{equateThresholds}{logical. Whether ordinal thresholds should be
   constrained equal across groups.}
}

\details{
For typical structural equation models the saturated model is the free-est possible model.  Not only all variances (and, when possible, all means) are estimated, but also all covariances.  In the case of ordinal data, the ordinal means are fixed to zero and the thresholds are estimated.  For binary variables, those variances are also constrained to one.  This is the free-est possible model that is identified.  The saturated model is used in calculating fit statistics such as the RMSEA, and Chi-squared fit indices.

The independence model, sometimes called the null model, is a model in which each variable is treated as being completely independent of every other variable.  As such, all the variances and, when possible, all means are estimated.  However, covariances are set to zero.  Ordinal variables are handled the same for the independence and saturated models.  The independence model is used, along with the saturated model, to create CFI and TLI fit indices.

The saturated and independence models could be used to create further fit indices.  However, OpenMx does not recommend using GFI, AGFI, NFI (aka Bentler-Bonett), or SRMR.  The page for \code{\link{mxSummary}} has information about why.

When the \link{mxFitFunctionMultigroup} fit function is used, \code{mxRefModels} creates the appropriate multi-group saturated and independence models.  Saturated and independence models are created separately for each group.  Each group has its own saturated and independence model.  The multi-group saturated model is a multi-group model where each group has its own saturated model, and similarly for the independence model.

When an MxModel has been run, some effort is made to make the reference models for only the variables used in the model.  For covariance data, all variables are modeled by default.  For raw data when the model has been run, only the modeled variables are used in the reference models.  This matches the behavior of \link{mxModel}.

In general, it is best practice to give \code{mxRefModels} a model that has already been run.

Multivariate normal models with all ordinal data and no missing values
can use the saturated multinomial distribution. This is much faster than
estimation of the saturated multivariate normal model. Use
\code{distribution='multinomial'} to avail this option.
}

\section{Warnings}{
One potentially important limitation of the \code{mxRefModels} function is for behavior-genetic models.  If variables 'x', 'y', and 'z' are measured on twins 1 and 2 creating the modeled variables 'x1', 'y1', 'z1', 'x2', 'y2', 'z2', then this function may not create the intended saturated or independence models.  In particular, the means of 'x1' and 'x2' are estimated separately.  Similarly, the covariance of 'x1' with 'y1' and 'x2' with 'y2' are allowed be be distinct: \eqn{cov(x1, y1) != cov{x2, y2}}.  Moreover, the cross-twin covariances are estimated: e.g. \eqn{cov(x1, y2) != 0}.

Another potential misuse of this function is for models with definition variables.  If definition variables are used, the saturated and independence model may not be correct because they do not account for the definition variables.

The are a few considerations specific to IFA models
(\link{mxExpectationBA81}).
The independence model preserves equality constraints among item
parameters from the original model.
The saturated model is a multinomial distribution with the proportions
equal to the proportions in your data. For example, if you have 2
dichotomous items then there are 4 possible response patterns: 00, 01, 10, 11.
A multinomial distribution for these 2 items is fully specified by 3
proportions or 3 parameters: a, b, c, \eqn{1.0-(a+b+c)}.
Hence, there is no need to optimize the saturated model.
When there is no missing data,
the deviance is immediately known as \eqn{-2 * sum(log proportions)}.
Typical Bayesian priors involve latent
factors (various densities on the pseudo-guessing lower bound, log norm
on loading, and uniqueness prior). These priors cannot be included
in the independence model because there are no latent factors.
Therefore, exercise caution when comparing the independence model
to a model that includes Bayesian priors.

\code{mxRefModels()} is not compatible with
\link[OpenMx:mxExpectationGREML]{GREML expectation}, as there is no
sensible general definition for a saturated GREML-type model.
}

\references{
The OpenMx User's guide can be found at http://openmx.ssri.psu.edu/documentation.
}

\examples{
require(OpenMx)
data(demoOneFactor)
manifests <- names(demoOneFactor)
latents <- c("G")
factorModel <- mxModel("OneFactor", type = "RAM",
    manifestVars = manifests, latentVars = latents,
    mxPath(from = latents, to=manifests, values = diag(var(demoOneFactor))*.2),
    mxPath(from = manifests, arrows = 2, values = diag(var(demoOneFactor))*.8),
    mxPath(from = latents, arrows = 2, free = FALSE, values = 1),
    mxData(cov(demoOneFactor), type = "cov", numObs = 500)
)
factorRun <- mxRun(factorModel)
factorSat <- mxRefModels(factorRun, run=TRUE)
summary(factorRun)
summary(factorRun, refModels=factorSat)

# A raw-data example where using mxRefModels adds fit indices
\donttest{
m1 <- mxModel("OneFactor", type = "RAM",
    manifestVars = manifests, latentVars = latents,
    mxPath(latents, to=manifests, values = diag(var(demoOneFactor))*.2),
    mxPath(manifests, arrows = 2, values = diag(var(demoOneFactor))*.8),
    mxPath(latents, arrows = 2, free = FALSE, values = 1),
    mxPath("one", to = latents, free = FALSE, values = 0),
    mxPath("one", to = manifests, values = 0),
    mxData(demoOneFactor, type = "raw")
)
m1 <- mxRun(m1)
summary(m1) # CFI, TLI, RMSEA missing
summary(m1, refModels=mxRefModels(m1, run = TRUE))
}

}

