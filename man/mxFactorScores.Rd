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

\name{mxFactorScores}
\alias{mxFactorScores}


\title{Estimate factor scores and standard errors}

\description{
This function creates the factor scores and their standard errors under different methods for an MxModel object that has an MxExpectationLISREL object.
}

\usage{
mxFactorScores(model, type=c('ML', 'WeightedML', 'Regression'))
}


\arguments{
   \item{model}{An MxModel object with an MxExpectationLISREL}
   \item{type}{The type of factor scores to compute}
}

\details{
This is a helper function to compute or estimate factor scores along with their standard errors.  The two maximum likelihood methods create a new model for each data row.  They then estimate the factor scores as free parameters in a model with a single data row.  For 'ML', the conditional likelihood of the data given the factor scores is optimized: \deqn{L(D|F)}.  For 'WeightedML', the joint likelihood of the data and the factor scores is optimized: \deqn{L(D, F) = L(D|F) L(F)}.  The WeightedML scores are akin to the empirical Bayes random effects estimates from mixed effects modeling.  They display the same kind of shrinkage as random effects estimates, and for the same reason: they account for the latent variable distribution in their estimation.  In many cases, especially for ordinal data or missing data, the weighted ML scores are to be preferred over alternatives (Estabrook & Neale, 2013).

For 'Regression', factor scores are computed based on a simple formula.  This formula is equivalent to the formula for the Kalman updated scores in a state space model with zero dynamics (Priestly & Subba Rao, 1975).  Thus, to compute the regression factor scores, the appropriate state space model is set-up and the mxKalmanScores function is used to produce the factor scores and their standard errors.
}

\value{
    An array with dimensions (Number of Rows of Data, Number of Latent Variables, 2).  The third dimension has the scores in the first slot and the standard errors in the second slot.
}

\seealso{
\link{mxKalmanScores}
}


\references{
Estabrook, R. & Neale, M. C.  (2013).  A Comparison of Factor Score Estimation Methods in the Presence of Missing Data: Reliability and an Application to Nicotine Dependence.  \emph{Multivariate Behavioral Research, 48,} 1-27.

Priestley, M. & Subba Rao, T.  (1975).  The estimation of factor scores and Kalman filtering for discrete parameter stationary processes.  \emph{International Journal of Control, 21}, 971-975.

The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{

# Create and estimate a factor model
require(OpenMx)
data(demoOneFactor)
manifests <- names(demoOneFactor)
latents <- c("G")
factorModel <- mxModel("OneFactor", 
                       type="LISREL",
                       manifestVars=list(exo=manifests), 
                       latentVars=list(exo=latents),
                       mxPath(from=latents, to=manifests),
                       mxPath(from=manifests, arrows=2),
                       mxPath(from=latents, arrows=2, free=FALSE, values=1.0),
                       mxPath(from='one', to=manifests),
                       mxData(observed=demoOneFactor[1:50,], type="raw"))
summary(factorRun <- mxRun(factorModel))

# Estimate factor scores for the model
r1 <- mxFactorScores(factorRun, 'Regression')


}
