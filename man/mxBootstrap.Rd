\name{mxBootstrap}
\alias{mxBootstrap}

\title{Repeatedly estimate model using resampling with replacement}

\description{
  Bootstrapping can help to assess the variability in parameter estimates.
  A new sample is draw from current data by unformly sampling the
  original data with replacement. Then the model is fit to these new
  data. This process is repeated many times. Quantiles of the estimates from all
  these replications can be used to assess the variability of parameters.
}

\usage{
mxBootstrap(model, replications=200, ...,
                        data=NULL, plan=NULL, verbose=0L,
                        parallel=TRUE, only=as.integer(NA),
			OK=c("OK", "OK/green"))
}

\arguments{
  \item{model}{The MxModel to be run.}
  \item{replications}{The number of resampling replications. If
    available, replications from prior mxBootstrap invocations will be reused.}
  \item{...}{Not used.  Forces remaining arguments to be specified by name.}
  \item{data}{A character vector of data or model names}
  \item{plan}{The compute plan used to optimize the model for each data
    set.}
  \item{verbose}{For levels greater than 0, enables runtime diagnostics}
  \item{parallel}{Whether to process the replications in parallel}
  \item{only}{When provided, only the given replication from a prior
    run of \code{mxBootstrap} will be performed. See details.}
  \item{OK}{The set of status code that are considered successful}
}

\details{
  By default, all datasets in the given model are resampled
  independently. If resampling is desired from only some of
  the datasets then the models containing them can be listed in the
  \sQuote{data} parameter.

  The \sQuote{plan} option defaults to using whatever compute plan is found in
  the given model (typically the output of
  \code{\link{omxDefaultComputePlan}}). Note, however, that the default
  plan often includes \code{\link{mxComputeNumericDeriv}}. This means
  that, after optimization, the Hessian will be approximated by finite
  differences. This procedure is of some value because it can be
  informative to check whether the Hessian is positive
  definite (see \code{\link{mxComputeHessianQuality}}).
  However, approximating the Hessian is often costly in terms
  of CPU time. For bootstrapping, the parameter estimates derived from
  the resampled data are typically of primary interest. Therefore, it is
  possible to specify a custom compute plan for bootstrap replications.
  For many models, \code{\link{mxComputeGradientDescent}} alone will suffice.

  On occasion, replications will fail. Sometimes it can be helpful to
  exactly reproduce a failed replication to attempt to pinpoint the
  cause of failure. The \sQuote{only} option facilitates this kind of
  investigation. In normal operation, mxBootstrap uses the regular R
  random number generator to generate a seed for each replication. This
  seed is used to seed an internal pseudorandom number generator
  (currently the Mersenne Twister algorithm). These
  per-replication seeds are stored as part of the bootstrap output. When
  \sQuote{only} is specified, the associated stored seed is use to seed the
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
  obtain per-parameter quantiles.
}

\seealso{
\code{\link{mxComputeBootstrap}}, \code{\link{mxSummary}}
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
boot <- mxBootstrap(boot)
summary(boot)

# examine replication 3
boot3 <- mxBootstrap(boot, only=3)

print(coef(boot3))
print(boot$compute$output$raw[3,names(coef(boot3))])
}
