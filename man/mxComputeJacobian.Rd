\name{mxComputeJacobian}
\alias{mxComputeJacobian}
\alias{MxComputeJacobian-class}
\title{Numerically estimate the Jacobian with respect to free parameters}
\usage{
mxComputeJacobian(freeSet=NA_character_, ..., of = "expectation", defvar.row=1L)
}
\arguments{
\item{freeSet}{names of matrices containing free variables}
\item{...}{Not used.  Forces remaining arguments to be specified by name.}
\item{of}{a character vector of expectations or algebra names}
\item{defvar.row}{A row index. Which row to load for definition variables.}
}
\description{
  When expectations are given, the Jacobian is taken with respect
  to the manifest model.
The manifest model excludes any latent variables or processes. For
RAM and LISREL models, the manifest model contains only the
manifest variables with free means, covariance, and thresholds.
Ordinal manifest variables are standardized.
}
\seealso{
  \link{omxManifestModelByParameterJacobian}, \link{mxGetExpected}
}
