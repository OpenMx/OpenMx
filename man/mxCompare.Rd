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

\name{mxCompare}
\alias{mxCompare}
\alias{mxCompareMatrix}

\title{Likelihood ratio test}

\description{
    Compare the fit of a model or set of models to a reference model or set of reference models.
}

\usage{
mxCompare(base, comparison, ..., all = FALSE,
  boot=FALSE, replications=400, previousRun=NULL, checkHess=FALSE)
mxCompareMatrix(models,
			    diag=c('minus2LL','ep','df','AIC'),
			    stat=c('p', 'diffLL','diffdf'), ...,
  boot=FALSE, replications=400, previousRun=NULL,
 checkHess=FALSE, wholeTable=FALSE)
}

\arguments{
   \item{base}{A MxModel object or list of MxModel objects.}
   \item{comparison}{A MxModel object or list of MxModel objects.}
   \item{models}{A MxModel object or list of MxModel objects.}
   \item{diag}{statistic used for diagonal entries}
   \item{stat}{statistic used for off-diagonal entries}
   \item{...}{Not used.}
   \item{all}{A boolean value on whether to compare all bases with all
     comparisons. Defaults to FALSE.}
   \item{boot}{Whether to use the bootstrap distribution to compute the
     p-value.}
   \item{replications}{How many replications to use to approximate the
     bootstrap distribution.}
   \item{previousRun}{Results to re-use from a previous bootstrap.}
  \item{checkHess}{Whether to approximate the Hessian in each replication}
  \item{wholeTable}{Return the whole table instead of a matrix shaped summary}
}

\details{
The mxCompare function is used to compare the fit of one or more
\link{MxMatrix} objects with output to one or more comparison
models. \code{mxCompareMatrix} compares all the models provided against
each other. Fit statistics for the comparison model or models are subtracted from the fit statistics for the base model or models. All models included in the \sQuote{base} argument are also listed without comparison (compared to a <NA> model) to present their raw fit statistics.

Model comparisons are made by subtracting the fit of the comparison
model from the fit of a base model. To make sure that the differences
between models are positive and yield p-values for likelihood ratio
tests, the model or models listed in the \sQuote{base} argument should
be more saturated (i.e., more estimated parameters and fewer degrees of
freedom) than models listed in the \sQuote{comparison} argument. For
\code{mxCompareMatrix} only the comparisons that make sense will be
included. If a comparison is made where the comparison model has a
higher minus 2 log likelihood (-2LL) than the base model, then the
difference in their -2LLs will be negative. P-values for likelihood
ratio tests will not be reported when either the -2LL or degrees of
freedom for the comparison are negative.

When multiple models are included in both the \sQuote{base} and \sQuote{comparison} arguments, then comparisons are made between the two lists of models based on the value of the \sQuote{all} argument. If \sQuote{all} is set to FALSE (default), then the first model in the \sQuote{base} list is compared to the first model in the \sQuote{comparison} list, second with second, and so on. If there are an unequal number of \sQuote{base} and \sQuote{comparison} models, then the shorter list of models is repeated to match the length of the longer list. For example, comparing base models \sQuote{B1} and \sQuote{B2} with comparison models \sQuote{C1}, \sQuote{C2} and \sQuote{C3} will yield three comparisons: \sQuote{B1} with \sQuote{C1}, \sQuote{B2} with \sQuote{C2}, and \sQuote{B1} with \sQuote{C3}. Each of those comparisons are prefaced by a comparison between the base model and a missing comparison model to present the fit of the base model.

If \sQuote{all} is set to TRUE, all possible comparisons between base and comparison models are made, and one entry is made for each base model. All comparisons involving the first model in \sQuote{base} are made first, followed by all comparisons with the second \sQuote{base} model, and so on. When there are multiple models in either the \sQuote{base} or \sQuote{comparison} arguments but not both, then the \sQuote{all} argument does not affect the set of comparisons made.

The following columns appear in the output:
\describe{
\item{base}{Name of the base model.}
\item{comparison}{Name of the comparison model. Is <NA> for the first }
\item{ep}{Estimated parameters of the comparison model.}
\item{minus2LL}{Minus 2*log-likelihood of the comparison model. If the comparison model is <NA>, then the minus 2*log-likelihood of the base model is given.}
\item{df}{Degrees in freedom of the comparison model. If the comparison model is <NA>, then the degrees of freedom of the base model is given.}
\item{AIC}{Akaike's Information Criterion for the comparison model. If the comparison model is <NA>, then the AIC of the base model is given.}
\item{diffLL}{Difference in minus 2*log-likelihoods of the base and comparison models. Will be positive when base model -2LL is higher than comparison model -2LL.}
\item{diffdf}{Difference in degrees of freedoms of the base and comparison models. Will be positive when base model DF is lower than comparison model DF (base model estimated parameters is higher than comparison model estimated parameters)}
\item{p}{P-value for likelihood ratio test based on diffLL and diffdf values.}}

The mxCompare function will give a p-value for any comparison in which
both \sQuote{diffLL} and \sQuote{diffdf} are non-negative. However, this
p-value is based on the assumptions of the likelihood ratio test,
specifically that the two models being compared are nested. The
likelihood ratio test and associated p-values are not valid when the
comparison model is not nested in the referenced base model. For a more
accurate p-value, the empirical bootstrap distribution can be computed
(\sQuote{boot=TRUE}). However, \sQuote{replications} must be set high
enough for an accurate approximation. The Monte Carlo SE of a proportion
for B replications is \eqn{\sqrt(p*(1-p)/B)}, but this will be zero if p
is zero, which is nonsense.  Note that a parametric-bootstrap p-value of
zero must be interpreted as \eqn{p < 1/B}, which, depending on \eqn{B}
and the desired Type I error rate, may not be "statistically significant."

When \sQuote{boot=TRUE}, the model has a default compute plan, and
\sQuote{checkHess} is kept at FALSE then the Hessian will not be
approximated or checked.  On the other hand, \sQuote{checkHess} is TRUE
then the Hessian will be approximated by finite differences. This
procedure is of some value because it can be informative to check
whether the Hessian is positive definite (see
\code{\link{mxComputeHessianQuality}}).  However, approximating the
Hessian is often costly in terms of CPU time. For bootstrapping, the
parameter estimates derived from the resampled data are typically of
primary interest.

Use options('digits' = N) to set the minimum number of significant digits to be printed in values. The mxCompare function does not directly accept a digits argument, and depends on the value of the 'digits' option.
}

\seealso{
\code{\link{mxPowerSearch}};
\code{\link{mxModel}};  \code{\link{options}} (use options('mxOptions') to see all the OpenMx-specific options) 
}

\examples{

data(demoOneFactor)
manifests <- names(demoOneFactor)
latents <- c("G1")
model1 <- mxModel(model="One Factor", type="RAM",
      manifestVars = manifests,
      latentVars = latents,
      mxPath(from = latents, to=manifests),
      mxPath(from = manifests, arrows = 2),
      mxPath(from = latents, arrows = 2, free = FALSE, values = 1.0),
      mxData(cov(demoOneFactor), type = "cov", numObs = 500)
)

fit1 <- mxRun(model1)

latents <- c("G1", "G2")
model2 <- mxModel(model="Two Factor", type="RAM",
      manifestVars = manifests,
      latentVars = latents,
      mxPath(from = latents[1], to=manifests[1:3]),
      mxPath(from = latents[2], to=manifests[4:5]),
      mxPath(from = manifests, arrows = 2),
      mxPath(from = latents, arrows = 2, free = FALSE, values = 1.0),
      mxData(cov(demoOneFactor), type = "cov", numObs=500)
)
fit2 <- mxRun(model2)

mxCompare(fit1, fit2)

# vary precision of the output
oldPrecision = as.numeric(options('digits')) 
options('digits' = 1)
mxCompare(fit1, fit2)
options('digits' = oldPrecision)
}
