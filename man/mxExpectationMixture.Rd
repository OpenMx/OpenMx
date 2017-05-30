\name{mxExpectationMixture}
\alias{MxExpectationMixture-class}
\alias{mxExpectationMixture}
\alias{print,MxExpectationMixture-method}
\alias{show,MxExpectationMixture-method}
\title{Mixture expectation}
\usage{
mxExpectationMixture(components, weights="weights",
		      ..., verbose=0L, scale=c('softmax','sum'))
}
\arguments{
  \item{components}{A character vector of model names.}

  \item{weights}{The name of the matrix or algebra column that specifies
  the component weights.}

\item{...}{Not used.  Forces remaining arguments to be specified by name.}

\item{verbose}{the level of runtime diagnostics}

\item{scale}{How the probabilities are rescaled. For 'softmax',
the coefficient-wise exponential is taken and then each column is
divided by its column sum. For 'sum', each column is divided by its
column sum.}
}
\description{
  Used in conjuction with \link{mxFitFunctionML}, this expectation
  can express a mixture model.
}
\details{
  Parameters are estimated in the given scale. To obtain the weights
  column vector, examine the expectation's \code{output} slot.

  An extension of this expectation to a Hidden Markov model
  is available with \link{mxExpectationHiddenMarkov}.
\link{mxGenerateData} is not implemented for this type of expectation.
}
\examples{
library(OpenMx)

set.seed(1)

trail <- c(rep(1,480), rep(2,520))
trailN <- sapply(trail, function(v) rnorm(1, mean=v))

classes <- list()

for (cl in 1:2) {
  classes[[cl]] <- mxModel(paste0("class", cl), type="RAM",
                           manifestVars=c("ob"),
                           mxPath("one", "ob", value=cl, free=FALSE),
                           mxPath("ob", arrows=2, value=1, free=FALSE),
                           mxFitFunctionML(vector=TRUE))
}

mix1 <- mxModel(
  "mix1", classes,
  mxData(data.frame(ob=trailN), "raw"),
  mxMatrix(values=1, nrow=1, ncol=2, free=c(FALSE,TRUE), name="weights"),
  mxExpectationMixture(paste0("class",1:2), scale="softmax"),
  mxFitFunctionML())

mix1Fit <- mxRun(mix1)
}
