\name{mxParametricBootstrap}
\alias{mxParametricBootstrap}

\title{Assess whether potential parameters should be freed using
  parametric bootstrap}

\description{
  Data is simulated from \sQuote{nullModel}. The parameters named by
  \sQuote{labels} are freed to obtain the alternative model. The
  alternative model is fit against each simulated data set.
}

\usage{
mxParametricBootstrap(nullModel, labels,
   alternative=c("two.sided", "greater", "less"),
   ..., alpha=0.05, correction=p.adjust.methods,
   previousRun=NULL, replications=400, checkHess=FALSE,
   signif.stars = getOption("show.signif.stars"))
}

\arguments{
  \item{nullModel}{The model specifying the null distribution}
  \item{labels}{A character vector of parameters to free to obtain
  the alternative model}
  \item{alternative}{a character string specifying the alternative hypothesis}
  \item{...}{Not used.  Forces remaining arguments to be specified by
    name.}
  \item{alpha}{The false positive rate}
  \item{correction}{How to adjust the p values for multiple tests.}
  \item{replications}{The number of resampling replications. If
    available, replications from prior invocation will be reused.}
  \item{previousRun}{Results to re-use from a previous bootstrap.}
  \item{checkHess}{Whether to approximate the Hessian in each
    replication}
  \item{signif.stars}{logical. If TRUE, \sQuote{significance stars} are
  printed for each coefficient.}
}

\details{
  When the model has a default compute plan and \sQuote{checkHess} is
  kept at FALSE then the Hessian will not be approximated or checked.
  On the other hand, \sQuote{checkHess} is TRUE then the Hessian will be
  approximated by finite differences. This procedure is of some value
  because it can be informative to check whether the Hessian is positive
  definite (see \code{\link{mxComputeHessianQuality}}).  However,
  approximating the Hessian is often costly in terms of CPU time. For
  bootstrapping, the parameter estimates derived from the resampled data
  are typically of primary interest.
}

\value{
A data frame is returned containing the test results.
Details of the bootstrap replications are stored in the
\sQuote{bootData} attribute on the data frame.
}

\examples{
library(OpenMx)

data(demoOneFactor)
manifests <- names(demoOneFactor)
latents <- c("G")

base <- mxModel(
  "OneFactorCov", type="RAM",
  manifestVars = manifests,
  latentVars = latents,
  mxPath(from=latents, to=manifests, values=0,  free=FALSE, labels=paste0('l',1:length(manifests))),
  mxPath(from=manifests, arrows=2, values=rlnorm(length(manifests)), lbound=.01),
  mxPath(from=latents, arrows=2, free=FALSE, values=1.0),
  mxPath(from = 'one', to = manifests, values=0, free=TRUE, labels=paste0('m',1:length(manifests))),
  mxData(demoOneFactor, type="raw"))

base <- mxRun(base)

# use more replications, 8 is for demonstration purpose only
mxParametricBootstrap(base, paste0('l', 1:length(manifests)),
                      "two.sided", replications=8)
}
