\name{mxPowerSearch}
\alias{mxPowerSearch}
\alias{mxPower}
\title{Power curve}
\description{
  Determine the power curve between true and false models
}
\usage{
mxPowerSearch(trueModel, falseModel, n=NULL, sig.level=0.05, ...,
                    probes=300L, previousRun=NULL,
                    gdFun=mxGenerateData,
                    method=c('empirical', 'ncp'),
                    grid=NULL,
                    statistic=c('LRT','AIC','BIC'),
		    OK=mxOption(trueModel, "Status OK"),
		checkHess=FALSE,
		silent=!interactive())

mxPower(trueModel, falseModel, n=NULL, sig.level=0.05, power=0.8, ...,
                    probes=300L, gdFun=mxGenerateData,
                    method=c('empirical', 'ncp'),
                    statistic=c('LRT','AIC','BIC'),
                    OK=mxOption(trueModel, "Status OK"), checkHess=FALSE,
                    silent=!interactive())
}
\arguments{
  \item{trueModel}{The true generating model for the data.}
  \item{falseModel}{The model representing the null hypothesis that we
    wish to reject.}
  \item{n}{Average group sample size. When NULL, finds the relationship
    between sample size and power. When given, finds the relationship
    between a parameter and power (at the given sample size).}
  \item{sig.level}{The false positive rate (a.k.a. type 1 error).}
  \item{power}{The power (a.k.a. 1 - type 2 error).}
  \item{...}{Not used.}
  \item{probes}{The number of probes to use when
    method=\sQuote{empirical}.}
  \item{previousRun}{Output from a prior run of \sQuote{mxPowerSearch} to extend.}
  \item{gdFun}{The function invoked to generate new data for each Monte
    Carlo probe.}
  \item{method}{To estimate power using the Monte Carlo method
    (\sQuote{empirical}) or the average non-centrality method
    (\sQuote{ncp}).}
  \item{grid}{A vector of locations at which to evaluate the power. If
    not provided, a reasonable default will be chosen.}
  \item{statistic}{Which test to use to compare models.}
  \item{OK}{The set of status code that are considered successful}
  \item{checkHess}{Whether to approximate the Hessian in each
    replication}
  \item{silent}{Whether to output a progress indicator.}
}
\details{

  Power is the chance of obtaining a significant difference when there
  is a significant difference (1 - false negative rate).
  The likelihood ratio test is used by default. There are two
  methods available to produce a power curve. The default,
  \code{method=empirical}, works for any model where the likelihood
  ratio test works. For example, definition variables and missing data
  are fine, but parameters estimated at upper or lower bounds will cause
  problems. The \code{method=empirical} can require a lot of time
  because the models need to be fit 100s of times. An alternate
  approach, \code{method=ncp}, is much more efficient and takes
  advantage of the fact that the non-null distribution of likelihood
  ratio test statistic is often \eqn{\chi^2(df_1 - df_0, N \lambda),
  chi^2(df1 - df0, N lambda)}. That is, the non-centrality parameter,
  \eqn{\lambda,lambda}, can be assumed, on average, to
  contribute equally per row. This permits essentially instant power curves
  without the burden of tedious simulation. However, definition
  variables, missing data, or unconventional models (e.g., mixture
  distributions) can violate this assumption. Therefore, we recommend
  verifying that the output from \code{method=empirical} roughly matches
  \code{method=ncp} on the model of interest before relying on
  \code{method=ncp}. Unlike \code{method=empirical}, \code{method=ncp}
  does not use the \code{gdFun} argument and works with models that
  only have summary statistics for data.

  When \code{method=ncp}, parameters of both \sQuote{trueModel} and
  \sQuote{falseModel} are assumed to be converged to
  their desired values.
  In contrast, when \code{method=empirical}, \sQuote{trueModel}
  need not be run or even contain data. On each replication,
  data are generated from \sQuote{trueModel} at the given parameter
  vector. Both \sQuote{trueModel} and \sQuote{falseModel} are
  fit against these data.

  When \code{statistic='LRT'} then the models must be nested and
  \code{sig.level} is used to determine whether the test is rejected.
  For \sQuote{AIC} and \sQuote{BIC}, the test is regarded as rejected
  when the \sQuote{trueModel} obtains a lower score than the
  \sQuote{falseModel}. In contrast to \code{statistic='LRT'}, there is
  no nesting requirement. For example, \sQuote{AIC} can be used to
  compare a \sQuote{trueModel} against its corresponding saturated
  model.

% discuss more than 1 parameter difference? (like an omnibus test)

\code{mxPower} operates in many modes.
When power is passed as \code{NULL} then power is calculated
and returned.
When power (as a scalar or vector) is given then sample or effect size
is (are) returned.
If you pass a list of models for \sQuote{falseModel}, each model
will be checked in turn and the results returned as a vector.
If you pass a vector of sample sizes, each sample size will
be checked in turn and the results returns as a vector.

In contrast to \code{mxPower}, \code{mxPowerSearch} attempts to model
the whole relationship between sample size or effect size and power.
  A naive Monte Carlo estimate of power examines a single candidate
  sample size at a time.  To obtain the whole curve, and simultaneously,
  to reduce the number of simulation probes, \code{mxPowerSearch}
  employs a binomial family generalized linear model with a logit link
  to predict the power curve across the sample sizes of interest
  (similar to Schoemann et al, 2014).

  The \code{mxPowerSearch} algorithm works in 3
  stages. Without loss of generality, our description will assume a
  sample size to power relationship, but a similar process is used for
  the parameter value to power relationship. In the first stage, a crude
  binary search is used to find the range reasonable values for N. This
  stage is complete once we have at least two rejections and at least two
  non-rejections. At this point, the binomial intercept and slope model
  is fit to these data. If the \emph{p}-value for the slope is less than
  0.25 then we jump to stage 3. Otherwise, we fit an intercept only
  binomial model. Our goal is to nail down the intercept (where
  power=0.5) because this is the easiest point to find and is a
  necessary prerequisite to estimate the variance (a.k.a. slope).
  Therefore, we probe at the median of previous probes stepping by 10\%
  in the direction of the model's predicted intercept. Eventually, after enough
  probes, we reach stage 3 where the \emph{p}-value for the slope is
  less than 0.25. At stage 3, our goal is to nail down the interesting
  part of the power curve. Therefore, we cycle serially through probes
  at 0, 1, and 2 logits from the intercept. This process is continued
  for the permitted number of \code{probes}.
  Occasionally, the \emph{p}-value for the slope in the stage 3 model
  grows larger than 0.25. In this case, we switch back to stage 2
  (intercept only) until the stage 3 model start working again.
  There is no other
  convergence criteria. Accuracy continues to improves until the probe
  limit is reached. After \code{mxPowerSearch} returns, additional probes can
  be run by passing the previous result back into the function as the
  \code{previousRun} argument.

  When \sQuote{n} is fixed then \code{mxPowerSearch} helps answer the
  question, \dQuote{how small of a true effect is likely to be detected
  at a particular sample size?} Only one parameter can be considered at
  a time.  The way the simulation works is that a candidate value for
  the parameter of interest is loaded into the \code{trueModel}, data
  are generated, then both the true and false model are fit to the data to
  obtain the difference in fit. The candidate parameter is initially set
  to halfway between the \code{trueModel} and \code{falseModel}.  The
  power curve will reflect the smallest distance, relative to the false
  model, required to have a good chance to reject the false model
  according to the chosen statistic.

  Note that the default \code{grid} is chosen to exhibit the interesting
  part of the power curve (from 0.25 to 0.97).
  Especially for \code{method=ncp}, this curve is practically
  identical for any pair of models (but located at a different range of
  sample sizes). You should select your own \code{grid} points if you
  wish to align power curves from more than one power analysis.

}

\value{
  \code{mxPower} returns a vector of sample sizes, powers, or effect sizes.

  \code{mxPowerSearch} returns a data.frame with one row for each \sQuote{grid} point. The
  first column is either the sample size \sQuote{N} or the parameter
  label. The second column is the power. The next two columns provide a
  +/-2 SE confidence interval for the power as estimated by the binomial
  logit linear model.
  When \code{method=empirical} then the \sQuote{probes} attribute
  contains a data.frame record of the activity of the power estimation
  process.

}

\references{

Schoemann, A. M., Miller, P., Pornprasertmanit, S. & Wu,
W. (2014). Using Monte Carlo simulations to determine power and sample
size for planned missing designs. \emph{International Journal of
  Behavioral Development, 38}(5), 471-479.

}

\seealso{
\code{\link{mxCompare}}, \code{\link{mxRefModels}}
}

\examples{
library(OpenMx)

data(demoOneFactor)
manifests <- names(demoOneFactor)
latents <- c("G")
factorModel <- mxModel("One Factor",
                       type="RAM",
                       manifestVars = manifests,
                       latentVars = latents,
                       mxPath(from=latents, to=manifests, values=0.8),
                       mxPath(from=manifests, arrows=2,values=1),
                       mxPath(from=latents, arrows=2,
                              free=FALSE, values=1.0),
                       mxPath(from="one", to=manifests),
                       mxData(demoOneFactor, type="raw"))
factorModelFit <- mxRun(factorModel)

indModel <- factorModelFit
indModel$A$values['x1','G'] <- 0.3
indModel$A$free['x1','G'] <- FALSE
indModel <- mxRun(indModel)

# What power do we have at different sample sizes
# to detect that the G to x1 factor loading is
# really 0.3 instead of 0.39?
mxPowerSearch(factorModelFit, indModel, method='ncp')

# If we want to conduct a study with 80% power to
# find that  the G to x1 factor loading is
# really 0.3 instead of 0.39, what sample size
# should we use?
mxPower(factorModelFit, indModel, method='ncp')
}
