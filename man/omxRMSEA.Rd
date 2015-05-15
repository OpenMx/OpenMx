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

\name{omxRMSEA}
\alias{omxRMSEA}


\title{Get the RMSEA with confidence intervals from model}

\description{
This function calculates the Root Mean Square Error of the Approximation (RMSEA) for a model and computes confidence intervals for that fit statistic.
}

\usage{
omxRMSEA(model, lower=.025, upper=.975, null=.05, ...)
}


\arguments{
   \item{model}{An MxModel object for which the RMSEA is desried}
   \item{lower}{The lower confidence bound for the confidence interval}
   \item{upper}{The upper confidence bound for the confidence interval}
   \item{null}{Value of RMSEA used to test for close fit}
   \item{...}{Further named arguments passed to summary}
}

\details{
To help users obtain fit statistics related to the RMSEA, this function confidence intervals and a test for close fit.  The user determines how close the fit is required to be by setting the \code{null} argument to the value desired for comparison.
}

\value{
    A named vector with elements lower, est.rmsea, upper, null, and `Prob(x <= null)`.
}


\references{
Browne, M. W. & Cudeck, R. (1992).  Alternative Ways of Assessing Model Fit.  \emph{Sociological Methods and Research, 21}, 230-258.
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
factorRun <- mxRun(factorModel)
factorSat <- mxRefModels(factorRun, run=TRUE)
summary(factorRun, refModels=factorSat)
# Gives RMSEA with 95% confidence interval

omxRMSEA(factorRun, .05, .95, refModels=factorSat)
# Gives RMSEA with 90% confidence interval
#  and probability of 'close enough' fit
}
