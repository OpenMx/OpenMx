\name{omxGetBootstrapReplications}
\alias{omxGetBootstrapReplications}
\alias{omxBootstrapCov}
\title{omxGetBootstrapReplications}
\usage{
omxGetBootstrapReplications(model)
omxBootstrapCov(model)
}
\arguments{
\item{model}{an MxModel object}
}
\value{
a matrix or covariance matrix of bootstrap parameter estimates
}
\description{
  Checks a variety of conditions to ensure that bootstrap replications are
  available and valid. Throws exception if things go wrong. Otherwise,
  replications are returned to the caller.
}
