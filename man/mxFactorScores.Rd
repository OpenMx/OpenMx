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

\name{mxFactorScores}
\alias{mxFactorScores}


\title{Estimate factor scores and standard errors}

\description{
This function creates the factor scores and their standard errors under different methods for an MxModel object that has either a RAM or LISREL expectation.
}

\usage{
mxFactorScores(model, type=c('ML', 'WeightedML', 'Regression'),
 minManifests=as.integer(NA))
}


\arguments{
   \item{model}{An MxModel object with either an MxExpectationLISREL or MxExpectationRAM}
   \item{type}{The type of factor scores to compute}
   \item{minManifests}{Set scores to NA when there are less than
   minManifests non-NA manifest variables}
}

\details{
This is a helper function to compute or estimate factor scores along with their standard errors.  The two maximum likelihood methods create a new model for each data row.  They then estimate the factor scores as free parameters in a model with a single data row.  For 'ML', the conditional likelihood of the data given the factor scores is optimized: \deqn{L(D|F)}.  For 'WeightedML', the joint likelihood of the data and the factor scores is optimized: \deqn{L(D, F) = L(D|F) L(F)}.  The WeightedML scores are akin to the empirical Bayes random effects estimates from mixed effects modeling.  They display the same kind of shrinkage as random effects estimates, and for the same reason: they account for the latent variable distribution in their estimation.

In many cases, especially for ordinal data or missing data, the weighted ML scores are to be preferred over alternatives (Estabrook & Neale, 2013). For example, when using ordinal data, a person whose observations are all in the highest ordinal category theoretically has an 'ML' factor score of positive infinity. A similar situation arises for a person whose observations are all in the lowest ordinal category: their 'ML' factor score is theoretically negative infinity. Weighted ML factor scores in these cases remain reasonable.

For type='Regression', with LISREL expectation, factor scores are computed based on a simple formula.  This formula is equivalent to the formula for the Kalman updated scores in a state space model with zero dynamics (Priestly & Subba Rao, 1975).  Thus, to compute the regression factor scores, the appropriate state space model is set-up and the \link{mxKalmanScores} function is used to produce the factor scores and their standard errors.  With RAM expectation, factor scores are predicted from the non-missing manifest variables for each row of the raw data, using a general linear prediction formula analytically equivalent to that used with LISREL expectation.  The standard errors for regression-predicted RAM factor scores are the square roots of the indeterminate variances of the latent variables, given the data row's missing-data pattern and the values of any relevant definition variables.  The RAM and LISREL methods for computing regression factor scores with their standard errors are analytically identical. They produce the same score and standard error estimates.

If you have missing data then you must specify \code{minManifests}.
This option will set scores to NA when there are too few items to make an accurate score estimate.
If you are using the scores as point estimates without considering the standard error then you should set minManifests as high as you can tolerate.
This will increase the amount of missing data but scores will be more accurate.
If you are carefully considering the standard errors of the scores then
you can set minManifests to 0.
When set to 0, all NA rows are scored to the prior distribution.

Note that for compatibility with \code{\link{factanal}}, type='regression' is also acceptable.
}

\value{
    An array with dimensions (Number of Rows of Data, Number of Latent Variables, 2).  The third dimension has the scores in the first slot and the standard errors in the second slot.  The rows are in the order of the \emph{unsorted} data.  Multigroup models are an exception, in that the returned value is instead a list of such arrays, containing one per group. 
}

\seealso{
\link{mxKalmanScores}
}


\references{
Estabrook, R. & Neale, M. C.  (2013).  A Comparison of Factor Score Estimation Methods in the Presence of Missing Data: Reliability and an Application to Nicotine Dependence.  \emph{Multivariate Behavioral Research, 48,} 1-27.

Priestley, M. & Subba Rao, T.  (1975).  The estimation of factor scores and Kalman filtering for discrete parameter stationary processes.  \emph{International Journal of Control, 21}, 971-975.

The OpenMx User's guide can be found at http://openmx.ssri.psu.edu/documentation.
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
                       mxData(observed=cov(demoOneFactor), type="cov", numObs=500,
                              means = colMeans(demoOneFactor)))
summary(factorRun <- mxRun(factorModel))

# Swap in raw data in place of summary data
factorRun <- mxModel(factorRun, mxData(observed=demoOneFactor[1:50,], type="raw"))

# Estimate factor scores for the model
r1 <- mxFactorScores(factorRun, 'Regression')


}
