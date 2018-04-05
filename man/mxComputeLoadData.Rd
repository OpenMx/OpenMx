\name{mxComputeLoadData}
\alias{mxComputeLoadData}
\alias{MxComputeLoadData-class}
\title{Load data from netcdf files directly into the backend}
\usage{
mxComputeLoadData(data, path)
}
\arguments{
\item{data}{a character vector of model names containing MxData objects}

\item{path}{a character vector of path templates}
}
\description{
Paths templates should be of the form \code{dataset-\%02d.nc}. This
compute step does nothing on first execution. Whatever data is initially
loaded is used. On second execution, it formats the path template with
the index count and loads these data from disk. The index counter starts
at 1 and is incremented every time data is loaded.
}
