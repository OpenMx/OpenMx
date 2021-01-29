%
%   Copyright 2007-2020 by the individuals mentioned in the source code history
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
    Compare the fit of one or more models to that of a reference (base) model or set of reference models.
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
   \item{all}{Boolean. Whether to compare all base models with all comparison models. Defaults to FALSE.}
   \item{boot}{Whether to use the bootstrap distribution to compute the p-value.}
   \item{replications}{How many replications to use to approximate the bootstrap distribution.}
   \item{previousRun}{Results to re-use from a previous bootstrap.}
  \item{checkHess}{Whether to approximate the Hessian in each replication}
  \item{wholeTable}{Return the whole table instead of a matrix shaped summary}
}

\details{
\code{mxCompare} is used to compare the fit of one or more \link{mxModel}s to one or more comparison models. \code{mxCompareMatrix} compares all the models provided against each other.

Model comparisons are made by subtracting the fit statistics for the comparison model from the fit statistics for the base model. Raw fit statistics of each \sQuote{base} model are also listed in the output table.

The fit statistics compared depend on the kinds of models compared.  Models fit with maximum likelihood are compared based on their minus two log likelihood values.  Under certain regularity conditions, the difference in minus two log likelihood values from nested models is chi-squared distributed and forms a likelihood ratio test statistic.  Models fit with weighted least squares are compared based on their Satorra-Bentler (2001) scaled difference chi-squared test statistics.  Under full weighted least squares, the Satorra-Bentler chi-squared value is equal to the difference in the model chi-squared values; however, for unweighted and diagonally weighted least squares, the two are no longer equal.  Satorra and Bentler (2001) showed that that their test statistic behaved well under a variety of conditions, including small sample sizes.  By contrast the much simpler difference in the chi-squared statistics only behaved well under large sample sizes (e.g., greater than or equal to 300 rows of data).

Specific to weighted least squares, researchers sometimes use mean-adjusted chi-squared statistics and mean-and-variance scaled chi-squared statistics.  Some programs call these WLSM and WLSMV statistics.  In some cases, it is fine to evaluate the total fit of a model using adjusted and scaled chi-squared statistics.  However, never, ever, ever, ..., ever take differences in mean-adjusted chi-squared statistics, and use them for nested model comparisons.  Similarly, never, ever, ever, ..., ever, ever take differences in mean-and-variance scaled chi-squared statistics, and use them for nested model comparisons.  The differences in these adjusted and scaled chi-squared statistics are not chi-squared distributed and do not form a valid basis for model comparison.  So, just don't do it.

Although not always checked by \code{mxCompare}, you should never compare models with different data sets or that use different variables from the same data set.  \code{mxCompare} might not stop you from doing this, so be thoughtful when comparing models.  Make sure your models are nested and use the same data.  Weighted least squares models are one case of comparing different data sets that requires particular care.  \emph{When comparing WLS models, make sure you are using the same exogenous covariates for all compared models.}  Because WLS is a multi-stage estimation approach, exogenous covariates residualize and change the data fitted in WLS.  Consequently, WLS models with different exogenous covariates actually have different data.  By contrast, maximum likelihood models with different exogenous covariates still use the same data and are valid to compare.

The \code{mxCompare} function makes an effort to only make valid comparisons.
If a comparison is made where the \code{comparison} model has a
higher minus 2 log likelihood (-2LL) than the \code{base} model, then the
difference in their -2LLs will be negative. P-values for likelihood
ratio tests will not be reported when either the -2LL or degrees of
freedom for the comparison are negative. 
To ensure that the differences between models are positive and yield p-values for likelihood ratio
tests, models listed in the \sQuote{base} argument must
be more saturated (i.e., more estimated parameters and fewer degrees of
freedom) than models listed in the \sQuote{comparison} argument. For
\code{mxCompareMatrix} only the comparisons that make sense will be included.

When multiple models are included in both the \sQuote{base} and \sQuote{comparison} arguments, then comparisons are made between the two lists of models based on the value of the \sQuote{all} argument. If \sQuote{all} is set to FALSE (default), then the first model in the \sQuote{base} list is compared to the first model in the \sQuote{comparison} list, second with second, and so on. If there are an unequal number of \sQuote{base} and \sQuote{comparison} models, then the shorter list of models is repeated to match the length of the longer list. For example, comparing base models \sQuote{B1} and \sQuote{B2} with comparison models \sQuote{C1}, \sQuote{C2} and \sQuote{C3} will yield three comparisons: \sQuote{B1} with \sQuote{C1}, \sQuote{B2} with \sQuote{C2}, and \sQuote{B1} with \sQuote{C3}. Each of those comparisons are prefaced by a comparison between the base model and a missing comparison model to present the fit of the base model.

If \sQuote{all} is set to TRUE, all possible comparisons between base and comparison models are made, and one entry is made for each base model. All comparisons involving the first model in \sQuote{base} are made first, followed by all comparisons with the second \sQuote{base} model, and so on. When there are multiple models in either the \sQuote{base} or \sQuote{comparison} arguments but not both, then the \sQuote{all} argument does not affect the set of comparisons made.

The following columns appear in the output for maximum likelihood comparisons:
\describe{
\item{base}{Name of the base model.}
\item{comparison}{Name of the comparison model. Is <NA> for the first }
\item{ep}{Estimated parameters of the comparison model.}
\item{minus2LL}{Minus 2*log-likelihood of the comparison model. If the comparison model is <NA>, then the minus 2*log-likelihood of the base model is given.}
\item{df}{Degrees in freedom of the comparison model. If the comparison model is <NA>, then the degrees of freedom of the base model is given.}
\item{AIC}{Akaike's An Information Criterion for the comparison model. If the comparison model is <NA>, then the AIC of the base model is given.}
\item{diffLL}{Difference in minus 2*log-likelihoods of the base and comparison models. Will be positive when base model -2LL is higher than comparison model -2LL.}
\item{diffdf}{Difference in degrees of freedoms of the base and comparison models. Will be positive when base model DF is lower than comparison model DF (base model estimated parameters is higher than comparison model estimated parameters)}
\item{p}{P-value for likelihood ratio test based on diffLL and diffdf values.}}

Weighted least squares reports a similar set of columns with four substitutions:
\describe{
\item{chisq}{Replaces the \code{minus2LL} column.  This is the comparison model's chi-squared statistic from Browne (1984, Equation 2.20a), accounting for some misspecification of the weight matrix.}
\item{AIC}{Although this has the same name as that in maximum likelihood, it is really a pseudo-AIC using the comparison model chi-squared and the number of estimated parameters.  It is the chi-squared value plus two times the number of free parameters.}
\item{SBchisq}{Replaces the \code{diffLL} column.  This is the Satorra-Bentler (2001, p. 511) scaled difference chi-squared statisic between the base model and the comparison model.  If your models use full weighted least squares, then this will be the same as the difference between the individual model chi-squared statistics.  However, for unweighted and diagonally weighted least square, the \code{SB chisq} will not be equal to the difference between the component model chi-squared statistics.}
\item{p}{p-value for the Satorra-Bentler chi-squared statistic.}}

In addition to the particular columns for maximum likelihood and weighted least squares, there are three general columns that are not printed but are accessible via the \code{$} and \code{[} extractors.
\describe{
\item{fit}{The individual model fit value: \code{m2ll} for maximum likelihood models, \code{chisq} for WLS models.}
\item{fitUnits}{The units of the fit function: \code{"-2LL"} for ML models, \code{"r'Wr"} for WLS models.}
\item{diffFit}{The difference in fit values between the base and comparison models: \code{diffLL} for ML models, \code{SBchisq} for WLS models.}}

\code{mxCompare} will give a p-value for any comparison in which
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

\emph{note}: The mxCompare function does not directly accept a digits argument, and depends
on the value of the 'digits' option. To set the minimum number of significant digits
printed, use options('digits' = N) (see example).
}

\value{
    Returns a new \code{MxCompare} object.  If you want something more like a table of results, use \code{as.data.frame()} on the returned \code{MxCompare} object.
}

\seealso{
\code{\link{mxPowerSearch}};
\code{\link{mxModel}};  \code{\link{options}} (use options('mxOptions') to see all the OpenMx-specific options) 
}

\examples{

data(demoOneFactor)
manifests <- names(demoOneFactor)
latents <- "G1"
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
model2 <- mxModel(model="One factor Rasch equated", type="RAM",
      manifestVars = manifests,
      latentVars = latents,
      mxPath(from = latents[1], to=manifests[1:5], labels='raschEquated'),

      mxPath(from = manifests, arrows = 2),
      mxPath(from = latents, arrows = 2, free = FALSE, values = 1.0),
      mxData(cov(demoOneFactor), type = "cov", numObs=500)
)
fit2 <- mxRun(model2)

mxCompare(fit1, fit2) # Rasch equated is significantly worse

# Vary precision (rounding) of the table 
oldPrecision = as.numeric(options('digits')) 
options('digits' = 1)
mxCompare(fit1, fit2)
options('digits' = oldPrecision)
}


