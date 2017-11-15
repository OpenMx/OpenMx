\name{mxPower}
\alias{mxPower}
\title{Power curve}
\description{
  Determine the power curve between true and false models using
 the likelihood ratio test
}
\usage{
mxPower(trueModel, falseModel, n=NULL, sig.level=0.05, ...,
                    probes=250L, prevOutput=NULL,
                    gdFun=mxGenerateData,
                    method=c('empirical', 'ncp'),
                    grid=NULL)
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
  \item{prevOutput}{Output from a prior run of \sQuote{mxPower} to extend.}
  \item{gdFun}{The function invoked to generate new data for each Monte
    Carlo probe.}
  \item{method}{To estimate power using the Monte Carlo method
    (\sQuote{empirical}) or the average non-centrality method
    (\sQuote{ncp}).}
  \item{grid}{A vector of locations at which to evaluate the power. If
    not provided, a reasonable default will be chosen.}
}
\details{
}

\value{
  Returns a data.frame with one row for each \sQuote{grid} point. The
  first column is either the sample size \sQuote{N} or the parameter
  label. The second column is the power. The next two columns provide a
  2 SE confidence interval for the power as estimated by the binomial
  logit linear model.

  When method=\sQuote{empirical} then the \sQuote{probes} attribute
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
