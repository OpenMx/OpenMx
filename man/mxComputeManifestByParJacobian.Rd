\name{mxComputeManifestByParJacobian}
\alias{mxComputeManifestByParJacobian}
\alias{MxComputeManifestByParJacobian-class}
\title{Manifest model by parameter Jacobian}
\usage{
mxComputeManifestByParJacobian(freeSet=NA_character_, ...,
 expectation = "expectation", defvar.row=1L)
}
\arguments{
\item{freeSet}{names of matrices containing free variables}
\item{...}{Not used.  Forces remaining arguments to be specified by name.}
\item{expectation}{a character vector of expectations to generate data
  for}
\item{defvar.row}{A row index. Which row to load for definition variables.}
}
\description{
The manifest model excludes any latent variables or processes. For
RAM and LISREL models, the manifest model contains only the
manifest variables with free means, covariance, and thresholds.
Ordinal manifest variables are standardized.
}
\seealso{
  \link{omxManifestModelByParameterJacobian}, \link{mxGetExpected}
}
