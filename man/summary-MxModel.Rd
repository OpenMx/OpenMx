\name{summary-MxModel}
\alias{summary,MxModel-method}

\title{Model Summary}

\description{
   This function returns summary statistics of a model. It is usually invoked after a model has been run through the optimizer.
}


\usage{
summary(model, ...)
}

\arguments{
   \item{model}{A MxModel object.}
   \item{...}{Any number of named arguments (see below).}
}

\details{

The following named arguments are supported by the summary method:
\itemize{
\item{SaturatedLikelihood}{Explicitly specify a saturated likelihood for testing (a numeric value)}
}

}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{

  model <- mxModel()
  modelOut <- mxRun(model)
  statistics <- mxSummary(modelOut)

}
