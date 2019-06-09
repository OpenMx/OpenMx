\name{mxBootstrap}
\alias{mxBootstrap}

\title{Repeatedly estimate model using resampling with replacement}

\description{
  Bootstrapping is used to quantify the variability of parameter estimates.
  A new sample is drawn from the model data (uniformly sampling the
  original data with replacement). The model is re-fitted to this new
  sample. This process is repeated many times. This yields a series of estimates
  from these replications which can be used to assess the variability of the parameters.
  
  \emph{note}: \code{mxBootstrap} only bootstraps free model parameters: 
  
  To bootstrap algebras, see \code{\link{mxBootstrapEval}}

  To report bootstrapped standardized paths in RAM models, \code{mxBootstrap} the model,
  and then run through \code{\link{mxBootstrapStdizeRAMpaths}}
  
}

\usage{
mxBootstrap(model, replications=200, ...,
                        data=NULL, plan=NULL, verbose=0L,
                        parallel=TRUE, only=as.integer(NA),
			OK=mxOption(model, "Status OK"), checkHess=FALSE)
}

\arguments{
  \item{model}{The MxModel to be run.}
  \item{replications}{The number of resampling replications. If
    available, replications from prior mxBootstrap invocations will be reused.}
  \item{...}{Not used.  Forces remaining arguments to be specified by name.}
  \item{data}{A character vector of data or model names}
  \item{plan}{Deprecated}
  \item{verbose}{For levels greater than 0, enables runtime diagnostics}
  \item{parallel}{Whether to process the replications in parallel (not yet implemented!)}
  \item{only}{When provided, only the given replication from a prior
    run of \code{mxBootstrap} will be performed. See details.}
  \item{OK}{The set of status code that are considered successful}
  \item{checkHess}{Whether to approximate the Hessian in each replication}
}

\details{
  By default, all datasets in the given model are resampled
  independently. If resampling is desired from only some of
  the datasets then the models containing them can be listed in the
  \sQuote{data} parameter.

  The \code{frequency} column in the \code{\link{mxData}} object is used
  represent a resampled dataset. When resampling, the original row
  proportions, as given by the original \code{frequency} column, are
  respected.

  When the model has a default compute plan and \sQuote{checkHess} is
  kept at FALSE then the Hessian will not be approximated or checked.
  On the other hand, \sQuote{checkHess} is TRUE then the Hessian will be
  approximated by finite differences. This procedure is of some value
  because it can be informative to check whether the Hessian is positive
  definite (see \code{\link{mxComputeHessianQuality}}).  However,
  approximating the Hessian is often costly in terms of CPU time. For
  bootstrapping, the parameter estimates derived from the resampled data
  are typically of primary interest.

  On occasion, replications will fail. Sometimes it can be helpful to
  exactly reproduce a failed replication to attempt to pinpoint the
  cause of failure. The \sQuote{only} option facilitates this kind of
  investigation. In normal operation, mxBootstrap uses the regular R
  random number generator to generate a seed for each replication. This
  seed is used to seed an internal pseudorandom number generator
  (currently the Mersenne Twister algorithm). These
  per-replication seeds are stored as part of the bootstrap output. When
  \sQuote{only} is specified, the associated stored seed is used to seed the
  internal random number generator so that identical weights can be
  regenerated.

  \sQuote{mxBootstrap} does not currently offer special support for nested or
  multilevel data. Rows are assumed model-wise independent.
}

\value{
  The given model is returned with
  the compute plan modified to consist of
  \code{mxComputeBootstrap}. Results of the bootstrap replications are
  stored inside the compute plan. \code{\link{mxSummary}} can be used to
  obtain per-parameter quantiles and standard errors.
}

\seealso{
  \code{\link{mxBootstrapEval}}, \code{\link{mxComputeBootstrap}},
  \code{\link{mxSummary}}, \code{\link{mxBootstrapStdizeRAMpaths}},
  \code{\link{as.statusCode}}
}

\examples{
library(OpenMx)

data(multiData1)

manifests <- c("x1", "x2", "y")

biRegModelRaw <- mxModel(
  "Regression of y on x1 and x2",
  type="RAM",
  manifestVars=manifests,
  mxPath(from=c("x1","x2"), to="y", 
         arrows=1, 
         free=TRUE, values=.2, labels=c("b1", "b2")),
  mxPath(from=manifests, 
         arrows=2, 
         free=TRUE, values=.8, 
         labels=c("VarX1", "VarX2", "VarE")),
  mxPath(from="x1", to="x2",
         arrows=2, 
         free=TRUE, values=.2, 
         labels=c("CovX1X2")),
  mxPath(from="one", to=manifests, 
         arrows=1, free=TRUE, values=.1, 
         labels=c("MeanX1", "MeanX2", "MeanY")),
  mxData(observed=multiData1, type="raw"))

biRegModelRawOut <- mxRun(biRegModelRaw)

boot <- mxBootstrap(biRegModelRawOut, 10)   # start with 10
summary(boot)

# Looks good, now do the rest
\donttest{boot <- mxBootstrap(boot)}
summary(boot)

# examine replication 3
boot3 <- mxBootstrap(boot, only=3)

print(coef(boot3))
print(boot$compute$output$raw[3,names(coef(boot3))])
}
