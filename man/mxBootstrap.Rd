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
mxBootstrap(model, replications=200, quantile=c(.25,.75), ...,
                        expectation=NULL, plan=NULL, verbose=0L,
                        parallel=TRUE, seed=42L, only=as.integer(NA),
			OK=c("OK", "OK/green"))
}

\arguments{
  \item{model}{The MxModel to be run.}
  \item{replications}{The number of resampling replications. Estimates from
    prior mxBootstrap invocations will be reused, if possible.}
  \item{quantile}{Quantiles used to summarize the parameter estimates}
  \item{...}{Not used.  Forces remaining arguments to be specified by name.}
  \item{expectation}{A character vector of expectation or model names}
  \item{plan}{The compute plan used to optimize the model for each data
    set.}
  \item{verbose}{For levels greater than 0, enables runtime diagnostics}
  \item{parallel}{Whether to process the replications in parallel}
  \item{seed}{The seed used for the internal random number generator}
  \item{only}{When provided, only the given replication will be
    performed. See details.}
  \item{OK}{The set of status code that are considered successful}
}

\details{
  If \sQuote{expectation} is 
}


\value{
Usually, \code{mxTryHard()} returns a post-\code{\link{mxRun}()} \code{\link{MxModel}} object.  Specifically, this will be the fitted model having the smallest fit-function value found by \code{mxTryHard()} during its attempts.  The start values used to obtain this fitted model are printed to console if \code{bestInitsOutput=TRUE}.  

If every attempt at running \code{model} fails, \code{mxTryHard()} returns an object of class 'try-error'.

\code{mxTryHard()} throws a warning if the returned \code{\link{MxModel}} object has a nonzero status code (unless nonzero status codes are considered acceptable per argument \code{greenOK} or \code{OKstatuscodes}).
}

\seealso{
\code{\link{mxRun}()}, \code{\link{mxComputeTryHard}}
}

\examples{ %Could use a better example...
library(OpenMx)

data(demoOneFactor)  # load the demoOneFactor dataframe

manifests <- names(demoOneFactor) # set the manifest to the 5 demo variables
latents   <- c("G")  # define 1 latent variable
model <- mxModel(model="One Factor", type="RAM",
    manifestVars = manifests,
    latentVars   = latents,
    mxPath(from=latents, to=manifests, labels=paste("b", 1:5, sep="")),
    mxPath(from=manifests, arrows=2, labels=paste("u", 1:5, sep="")),
    mxPath(from=latents  , arrows=2, free=FALSE, values=1.0),
    mxData(cov(demoOneFactor), type="cov", numObs=500)
)
model <- mxTryHard(model) # Run the model, returning the result into model
summary(model) # Show summary of the fitted model

}
