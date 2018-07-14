\name{mxComputeLoadData}
\alias{mxComputeLoadData}
\alias{write.cbor}
\alias{MxComputeLoadData-class}
\title{Load data from netcdf files directly into the backend}
\usage{
mxComputeLoadData(dest, path, ..., originalDataIsIndexOne=FALSE)
write.cbor(mxd, file)
}
\arguments{
\item{dest}{a character vector of model names containing MxData objects}

\item{path}{a character vector of path templates}

\item{...}{Not used.  Forces remaining arguments to be specified by name.}

\item{originalDataIsIndexOne}{logical. Whether to use the initial data
for index 1}

\item{mxd}{an mxData object}
\item{file}{a character string naming a file}
}
\description{
  Use \code{write.cbor} to store mxData objects containing datasets
  to files in cbor format.

  Paths templates should be of the form \code{dataset-\%02d.cbor}.
  The innermost \link{mxComputeLoop} index is substituted into this
  template using a printf-style format.
  During normal operation, the path template is formatted
  with the index count and data are loaded from disk.

  If \code{originalDataIsIndexOne} is TRUE then this
  compute step does nothing when the loop index is 1.
  The purpose of \code{originalDataIsIndexOne} is to
  permit usage of the dataset that was initially
  included with the model.
}
