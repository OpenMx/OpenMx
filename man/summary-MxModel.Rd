\name{summary-MxModel}
\alias{summary,MxModel-method}
\alias{summary}

\title{Model Summary}

\description{
   This function returns summary statistics of a model. It is usually invoked after a model has been run through the optimizer.
}


\usage{
summary(object, ...)
}

\arguments{
   \item{object}{A MxModel object.}
   \item{...}{Any number of named arguments (see below).}
}

\details{

The following named arguments are supported by the summary method:
\describe{
\item{SaturatedLikelihood}{Explicitly specify a saturated likelihood for testing (a numeric value)}
\item{indep}{Set to FALSE to ignore independent submodels in summary.}
}

}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{

  model <- mxModel()
  modelOut <- mxRun(model)
  # compute a summary and store in variable "statistics"
  statistics <- summary(modelOut)

  # specify a saturated likelihood for testing
  summary(modelOut, SaturatedLikelihood=300)

}
