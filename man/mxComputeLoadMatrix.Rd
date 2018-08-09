\name{mxComputeLoadMatrix}
\alias{mxComputeLoadMatrix}
\alias{MxComputeLoadMatrix-class}
\title{Load data from CSV files directly into the backend}
\usage{
mxComputeLoadMatrix(dest, path, ..., originalDataIsIndexOne=FALSE,
 row.names=FALSE, col.names=FALSE)
}
\arguments{
\item{dest}{a character vector of matrix names}

\item{path}{a character vector of paths}

\item{...}{Not used.  Forces remaining arguments to be specified by name.}

\item{originalDataIsIndexOne}{logical. Whether to use the initial data
for index 1}
\item{row.names}{logical. Whether row names are present}
\item{col.names}{logical. Whether column names are present}
}
\description{
  The CSV file must be formatted in a specific way.
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

Cov <- rWishart(4, 20, toeplitz(c(2,1)/20))
write.table(t(apply(Cov, 3, vech)),
            file="cov.csv", col.names=FALSE, row.names=FALSE)
Mean <- matrix(rnorm(8),4,2)
write.table(Mean, file="mean.csv", col.names=FALSE, row.names=FALSE)

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
                        c('mean.csv', 'cov.csv')),
    mxComputeOnce('fitfunction', 'fit'),
    mxComputeCheckpoint(path="loadMatrix.csv")
  ), i=1:4))

m1 <- mxRun(m1)
}
