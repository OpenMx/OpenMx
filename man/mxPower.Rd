\name{mxPower}
\alias{mxPower}
\title{Power curve}
\description{
  Determine the power curve between true and false models using
 the likelihood ratio test
}
\usage{
mxPower(trueModel, falseModel, n=NULL, sig.level=0.05, ...,
                    probes=250L, previousRun=NULL,
                    gdFun=mxGenerateData,
                    method=c('empirical', 'ncp'),
                    grid=NULL,
		    OK=mxOption(trueModel, "Status OK"), checkHess=FALSE)
}
\arguments{
  \item{trueModel}{The true generating model for the data.}
  \item{falseModel}{The model representing the null hypothesis that we
    wish to reject.}
  \item{n}{Average group sample size. When NULL, finds the relationship
    between sample size and power. When given, finds the relationship
    between a parameter and power (at the given sample size).}
  \item{sig.level}{The false positive rate (a.k.a. type 1 error).}
  \item{...}{Not used.}
  \item{probes}{The number of probes to use when
    method=\sQuote{empirical}.}
  \item{previousRun}{Output from a prior run of \sQuote{mxPower} to extend.}
  \item{gdFun}{The function invoked to generate new data for each Monte
    Carlo probe.}
  \item{method}{To estimate power using the Monte Carlo method
    (\sQuote{empirical}) or the average non-centrality method
    (\sQuote{ncp}).}
  \item{grid}{A vector of locations at which to evaluate the power. If
    not provided, a reasonable default will be chosen.}
  \item{OK}{The set of status code that are considered successful}
  \item{checkHess}{Whether to approximate the Hessian in each replication}
}
\details{

  To reduce the number of simulation probes,
  \code{method=\sQuote{empirical}} ultimately employs a binomial family
  generalized linear model with a logit link to predict the power curve
  across the sample sizes of interest. The algorithm works in 3
  stages. Without loss of generality, our description will assume a
  sample size to power relationship, but a similar process is used for
  the parameter value to power relationship. In the first stage, a crude
  binary search is used to find the range reasonable values for N. This
  stage is complete once we have at least one rejection and at least one
  non-rejection. At this point, the binomial intercept and slope model
  is fit to these data. If the \emph{p}-value for the slope is less than
  0.25 then we jump to stage 3. Otherwise, we fit an intercept only
  binomial model. Our goal is to nail down the intercept (where
  power=0.5) because this is the easiest point to find and is a
  necessary prerequisite to estimate the variance (or slope).
  Therefore, we probe at the median of previous probes stepping by 10\%
  in the direction of the model's intercept. Eventually, after enough
  probes, we reach stage 3 where the \emph{p}-value for the slope is
  less than 0.25. At stage 3, our goal is to nail down the interesting
  part of the power curve. Therefore, we cycle serially through probes
  at 0, 1, and 2 logits from the intercept. This process is continued
  for the permitted number of \code{probes}. There is no other
  convergence criteria. Accuracy continues to improves until the probe
  limit is reached. After \code{mxPower} returns, additional probes can
  be run by passing the previous result back into the function as the
  \code{previousRun} argument.

  Note that the default \code{grid} is chosen to exhibit the interesting
  part of the power curve (from 0.25 to 0.97 for \code{sig.level} of
  0.05). Especially for \code{method=\sQuote{ncp}}, this curve is practically
  identical for any pair of models (but located at a different range of
  sample sizes). You should select your own \code{grid} points if you
  wish to compare power curves from more than one power analyses.

}

\value{
  Returns a data.frame with one row for each \sQuote{grid} point. The
  first column is either the sample size \sQuote{N} or the parameter
  label. The second column is the power. The next two columns provide a
  2 SE confidence interval for the power as estimated by the binomial
  logit linear model.

  When \code{method=\sQuote{empirical}} then the \sQuote{probes} attribute
  contains a data.frame record of the activity of the power estimation
  process.

}

\references{

Schoemann, A. M., Miller, P., Pornprasertmanit, S. & Wu,
W. (2014). Using Monte Carlo simulations to determine power and sample
size for planned missing designs. \emph{International Journal of
  Behavioral Development, 38}(5), 471-479.

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

mxPower(factorModelFit, indModel, method='ncp')
}
