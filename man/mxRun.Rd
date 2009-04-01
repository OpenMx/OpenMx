\name{mxRun}
\alias{mxRun}

\title{Send a Model to the Optimizer}

\description{
   This function begins optimization on the top-level model.
}

\usage{
mxRun(model)
}

\arguments{
    \item{model}{An \link{MxModel} object}
}

\value{
   Return a model with free parameters updated to their final values.
   The return value contains an "output" slot with the results of 
   optimization.
}

\examples{

   model <- mxModel()

   model <- mxRun(model)

}
