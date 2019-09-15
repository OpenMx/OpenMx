\name{mxComputeLoadMatrix}
\alias{mxComputeLoadMatrix}
\alias{MxComputeLoadMatrix-class}
\title{Load data from CSV files directly into the backend}
\usage{
mxComputeLoadMatrix(dest, method=c('csv','data.frame'), ...,
 path=NULL, originalDataIsIndexOne=FALSE,
 row.names=FALSE, col.names=FALSE, observed=NULL)
}
\arguments{
\item{dest}{a character vector of matrix names}

\item{method}{name of the conduit used to load the data.}

\item{...}{Not used.  Forces remaining arguments to be specified by name.}

\item{path}{a character vector of paths}

\item{originalDataIsIndexOne}{logical. Whether to use the initial data
for index 1}
\item{row.names}{logical. Whether row names are present}
\item{col.names}{logical. Whether column names are present}
\item{observed}{data frame. The reservoir of data for \code{method='data.frame'}}
}
\description{
  THIS INTERFACE IS EXPERIMENTAL AND SUBJECT TO CHANGE.
  
  For method='csv', the file must be formatted in a specific way.
  The number of columns must match the number of entries available
  in the mxMatrix. Matrix types (e.g., symmetric or diagonal) are
  respected (see \link{mxMatrix}). For example, a \emph{Full} 2x2
  matrix will require 4 entries, but a diagonal matrix of the same size
  will only require 2 entries.
  CSV data must be stored space separated and without row or column
  names.
  The destination \code{mxMatrix} can have free parameters, but cannot
  have square bracket populated entries.

  If \code{originalDataIsIndexOne} is TRUE then this
  compute step does nothing when the loop index is 1.
  The purpose of \code{originalDataIsIndexOne} is to
  permit usage of the dataset that was initially
  included with the model.
}
\seealso{
  \link{mxComputeLoadData}, \link{mxComputeCheckpoint}
}
\examples{
library(OpenMx)

dir <-tempdir()  # safe place to create files

Cov <- rWishart(4, 20, toeplitz(c(2,1)/20))
write.table(t(apply(Cov, 3, vech)),
            file=file.path(dir, "cov.csv"),
            col.names=FALSE, row.names=FALSE)
Mean <- matrix(rnorm(8),4,2)
write.table(Mean, file=file.path(dir, "mean.csv"),
            col.names=FALSE, row.names=FALSE)

m1 <- mxModel(
  "test1",
  mxMatrix("Full", 1,2, values=0,       name="mean"),
  mxMatrix("Symm", 2,2, values=diag(2), name="cov"),
  mxMatrix("Full", 1,2, values=-1,      name="lbound"),
  mxMatrix("Full", 1,2, values=1,       name="ubound"),
  mxAlgebra(omxMnor(cov,mean,lbound,ubound), name="area"),
  mxFitFunctionAlgebra("area"),
  mxComputeLoop(list(
    mxComputeLoadMatrix(c('mean', 'cov'),
                        path=file.path(dir, c('mean.csv', 'cov.csv'))),
    mxComputeOnce('fitfunction', 'fit'),
    mxComputeCheckpoint(path=file.path(dir, "loadMatrix.csv"))
  ), i=1:4))

m1 <- mxRun(m1)
}
