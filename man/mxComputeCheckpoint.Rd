\name{mxComputeCheckpoint}
\alias{mxComputeCheckpoint}
\alias{MxComputeCheckpoint-class}
\title{Log parameters and state to disk or memory}
\usage{
mxComputeCheckpoint(what, ..., path=NULL, toReturn=TRUE,
				parameters=TRUE, loopIndices=TRUE, fit=TRUE, counters=TRUE)
}
\arguments{
\item{what}{a character vector of algebra names to include in each checkpoint}

\item{...}{Not used.  Forces remaining arguments to be specified by name.}

\item{path}{a character vector of where to write the checkpoint file}

\item{toReturn}{logical. Whether to store the checkpoint in memory and
  return it after the model is run}

\item{parameters}{logical. Whether to include the parameter vector}
\item{loopIndices}{logical. Whether to include the loop indices}
\item{fit}{logical. Whether to include the fit value}
\item{counters}{logical. Whether to include counters (number of
  evaluations and iterations)}
}
\description{
Captures the current state of the backend. When \code{path} is set, the
state is written to disk in a single row. When \code{toReturn} is set,
the state is recorded in memory and returned after \code{mxRun}.
}
