\name{mxComputeCheckpoint}
\alias{mxComputeCheckpoint}
\alias{MxComputeCheckpoint-class}
\title{Log parameters and state to disk or memory}
\usage{
mxComputeCheckpoint(what=NULL, ..., path=NULL, append=FALSE, header=TRUE, toReturn=FALSE,
				parameters=TRUE, loopIndices=TRUE,
  fit=TRUE, counters=TRUE, status=TRUE, standardErrors=FALSE,
  gradient=FALSE, vcov=FALSE)
}
\arguments{
\item{what}{a character vector of algebra names to include in each checkpoint}

\item{...}{Not used.  Forces remaining arguments to be specified by name.}

\item{path}{a character vector of where to write the checkpoint file}

\item{append}{if FALSE, truncates the checkpoint file upon open. If
  TRUE, existing data is preserved and checkpoints are appended.}

\item{header}{whether to write the header that describes the content of
  each column}

\item{toReturn}{logical. Whether to store the checkpoint in memory and
  return it after the model is run}

\item{parameters}{logical. Whether to include the parameter vector}
\item{loopIndices}{logical. Whether to include the loop indices}
\item{fit}{logical. Whether to include the fit value}
\item{counters}{logical. Whether to include counters (number of
  evaluations and iterations)}
\item{status}{logical. Whether to include the status code}
\item{standardErrors}{logical. Whether to include the standard errors}
\item{gradient}{logical. Whether to include the gradients}
\item{vcov}{logical. Whether to include the vcov in half-vectorized order}
}
\description{
Captures the current state of the backend. When \code{path} is set, the
state is written to disk in a single row. When \code{toReturn} is set,
the state is recorded in memory and returned after \code{mxRun}.
}
\examples{
library(OpenMx)

m1 <- mxModel(
  "poly22", # Eqn 22 from Tsallis & Stariolo (1996)
  mxMatrix(type='Full', values=runif(4, min=-1e6, max=1e6),
           ncol=1, nrow=4, free=TRUE, name='x'),
  mxAlgebra(sum((x*x-8)^2) + 5*sum(x) + 57.3276, name="fit"),
  mxFitFunctionAlgebra('fit'))

plan <- mxComputeLoop(list(
  mxComputeSetOriginalStarts(),
    mxComputeSimAnnealing(method="tsallis1996",
                          control=list(tempEnd=1)),
    mxComputeCheckpoint(path = "result.log")),
  i=1:4)

m1 <- mxRun(mxModel(m1, plan)) # see the file 'result.log'
}
