\name{mxComputeBootstrap}
\alias{mxComputeBootstrap}
\alias{MxComputeBootstrap-class}

\title{Repeatedly estimate model using resampling with replacement}

\description{This is a low-level compute plan object to perform
  resampling with replacement. }

\usage{
mxComputeBootstrap(data, plan, replications=200, ...,
                        verbose=0L, parallel=TRUE, freeSet=NA_character_,
			OK=c("OK", "OK/green"), only=NA_integer_)
}

\arguments{
  \item{data}{A vector of dataset or model names.}
  \item{plan}{The compute plan used to optimize the model for each data
    set.}
  \item{replications}{The number of resampling replications. If
    available, replications from prior mxBootstrap invocations will be reused.}
  \item{...}{Not used.  Forces remaining arguments to be specified by name.}
  \item{verbose}{For levels greater than 0, enables runtime diagnostics}
  \item{parallel}{Whether to process the replications in parallel}
  \item{freeSet}{names of matrices containing free variables}
  \item{OK}{The set of status code that are considered successful}
  \item{only}{When provided, only the given replication from a prior
    run of \code{mxBootstrap} will be performed. See details.}
}

\details{
  The \sQuote{only} option facilitates investigation of a single
  replication attempt. 
}


\value{
  Output is stored in the compute object's \code{output}
  slot. Specifically, \code{model$compute$output$raw} contains a data
  frame with parameters in columns and replications in rows. In addition
  to parameters, the \code{seed}, \code{fit}, and \code{statusCode}
  of the replication is also included.

  When \sQuote{only} is set to a particular replications, the weight
  vectors (one per dataset) are also returned in the compute object's
  \code{output} slot. \code{model$compute$output$weight} is a character
  vector (by dataset name) of numeric vectors (the weights). These
  weights can be used to recreate a model identical to the model used
  in the given replication.
}

\seealso{
\code{\link{mxBootstrap}}
}
