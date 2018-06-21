\name{mxExpectationMixture}
\alias{MxExpectationMixture-class}
\alias{mxExpectationMixture}
\alias{print,MxExpectationMixture-method}
\alias{show,MxExpectationMixture-method}
\title{Mixture expectation}
\usage{
mxExpectationMixture(components, weights="weights",
		      ..., verbose=0L, scale=c('softmax', 'sum', 'none'))
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
column sum. For 'none', no scaling is done.}
}
\description{
  Used in conjuction with \link{mxFitFunctionML}, this expectation
  can express a mixture model.
}
\details{
  The mixture probabilities given in \code{weights} must sum to one.  As such for \eqn{K} mixture components, only \eqn{K-1} of the elements of \code{weights} can be estimated. The mixture probabilities in \code{weights} should be a column vector (i.e., a \eqn{K} by 1 matrix, or algebra with a \eqn{K} by 1 result).
  
  For ease of use the raw free parameters of weights can be rescaled by OpenMx according to the \code{scale} argument.  When \code{scale} is set to "softmax" the softmax function is applied to the weights.  The softmax function is also sometimes called multinomial logistic regression.  Softmax exponentiates each element in a vector and then divides each element by the sum of the exponentiated elements.  In equation form the softmax function is
  
  \deqn{ softmax(x_i) = \frac{e^{x_i}}{\sum_{k=1}^{K} } e^{x_k} }{
    softmax(x_i) = exp(x_i) / sum(exp(x_i))}
  
  When using the softmax scaling no free parameter bounds or constraints are needed.  However, for model identification, one element of the weights vector must be fixed.  If the softmax scaling is used, then the usual choise for the fixed parameter value is zero.  The latent class or mixture component that has its raw weight set to zero becomes the comparison against which other probabilities are evaluated.
  
  When \code{scale} is set to "sum" then each element of the weights matrix is internally divided by its sum.  When using the sum scaling, the same model identification requirements are present.  In particular, one element of the weights must be fixed.  The typical value to fix this value at for sum scaling is one.  Additionally when using sum scaling, all free parameters in the weights must have lower bounds of zero.  In equation form the sum scaling does the following:
  
  \deqn{ sumscale(x_i) = \frac{x_i}{\sum_{k=1}^{K} } x_k }{
    sumscale(x_i) = x_i / sum(x_i)}
  
  When \code{scale} is set to "none" then no re-scaling is done. The weights are left "as is".  This can be dangerous and is not recommended for novice users.  However, some advanced users may find no scaling to be advantageous for certain applications (e.g., they are providing their own scaling), and thus it is provided as an option.
  
  Parameters are estimated in the given scale. To obtain the weights
  column vector, examine the expectation's \code{output} slot with for example \code{yourModel$expectation$output}

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
