\name{mxExpectationHiddenMarkov}
\alias{MxExpectationHiddenMarkov-class}
\alias{mxExpectationHiddenMarkov}
\alias{print,MxExpectationHiddenMarkov-method}
\alias{show,MxExpectationHiddenMarkov-method}
\title{Hidden Markov expectation}
\usage{
mxExpectationHiddenMarkov(components, initial="initial", transition=NULL,
		      ..., verbose=0L, scale=c('softmax', 'sum', 'none'))
}
\arguments{
  \item{components}{A character vector of model names.}

  \item{initial}{The name of the matrix or algebra column that specifies
  the initial probabilities.}

  \item{transition}{The name of the matrix or algebra that specifies
  the left stochastic transition probabilities.}

\item{...}{Not used.  Forces remaining arguments to be specified by name.}

\item{verbose}{the level of runtime diagnostics}

\item{scale}{How the probabilities are rescaled. For 'softmax',
the coefficient-wise exponential is taken and then each column is
divided by its column sum. For 'sum', each column is divided by its
column sum. For 'none', no scaling is done.}
}
\description{
  Used in conjunction with \link{mxFitFunctionML}, this expectation
  can express a mixture model (with the transition matrix omitted)
  or a Hidden Markov model.
}
\details{
  The initial probabilities given in \code{initial} must sum to one.  So too must the \emph{columns} of the transition matrix given in \code{transition}.  The transitions go from a column to a row, similar to how regression effects in the RAM structural equation models go from the column variable to the row variable.  This means \code{transition} is a left stochastic matrix.
  
  For ease of use the raw free parameters of these matrices are rescaled by OpenMx according to the \code{scale} argument.  When \code{scale} is set to "softmax" the softmax function is applied to the initial probabilities and the columns of the transition matrix.  The softmax function is also sometimes called multinomial logistic regression.  Softmax exponentiates each element in a vector and then divides each element by the sum of the exponentiated elements.  In equation form the softmax function is
  
  \deqn{ softmax(x_i) = \frac{e^{x_i}}{\sum_{k=1}^{K} } e^{x_k} }{
    softmax(x_i) = exp(x_i) / sum(exp(x_i))}
  
  When using the softmax scaling no free parameter bounds or constraints are needed.  However, for model identification, one element of the initial probabilities vector must be fixed.  If the softmax scaling is used, then the usual choice for the fixed parameter value is zero.  The regime (or latent class or mixture component) that has its initial probability set to zero becomes the comparison against which other probabilities are evaluated.  Likewise for model identification, one element in each column of the transition matrix must be fixed.  When the softmax scaling is used, the typical choice is to fix one element in each column to zero.  Generally, one row of the transition matrix is fixed to zero, or the diagonal elements of the transition matrix are fixed to zero.
  
  When \code{scale} is set to "sum" then each element of the initial probabilities and each column of the transition matrix is internally divided by its sum.  When using the sum scaling, the same model identification requirements are present.  In particular, one element of the initial probabilities must be fixed and one element in each column of the transition matrix must be fixed.  The typical value to fix these values at for sum scaling is one.  Additionally when using sum scaling, all free parameters in the initial and transition probabilities must have lower bounds of zero.  In equation form the sum scaling does the following:
  
  \deqn{ sumscale(x_i) = \frac{x_i}{\sum_{k=1}^{K} } x_k }{
    sumscale(x_i) = x_i / sum(x_i)}
  
  When \code{scale} is set to "none" then no re-scaling is done. The parameters of \code{initial} and \code{transition} are left "as is".  This can be dangerous and is not recommended for novice users.  It might produce nonsensical results particularly for hidden Markov models.  However, some advanced users may find no scaling to be advantageous for certain applications (e.g., they are providing their own scaling), and thus it is provided as an option.
  
  Parameters are estimated in the given scale. To obtain the initial
  column vector and left stochastic transition matrix in probability
  units then examine the expectation's \code{output} slot with for example \code{yourModel$expectation$output}

  Definition variables can be used to assign a separate set of mixture
  probabilities to each row of data.
  Definition variables can be used in the initial column vector or in
  the transition matrix, but not in both at the same time.

  Note that, when the transition matrix is omitted,
  this expectation is the same as \link{mxExpectationMixture}.
  \link{mxGenerateData} is not implemented for this type of expectation.
}
\examples{
library(OpenMx)

start_prob <- c(.2,.4,.4)
transition_prob <- matrix(c(.8, .1, .1,
			.3, .6, .1,
			.1, .3, .6), 3, 3)
noise <- .5

  # simulate a trajectory
  state <- sample.int(3, 1, prob=transition_prob \%*\% start_prob)
  trail <- c(state)

  for (rep in 1:500) {
    state <- sample.int(3, 1, prob=transition_prob[,state])
    trail <- c(trail, state)
  }
  
  # add noise
  trailN <- sapply(trail, function(v) rnorm(1, mean=v, sd=sqrt(noise)))

  classes <- list()
  
  for (cl in 1:3) {
    classes[[cl]] <- mxModel(paste0("cl", cl), type="RAM",
                             manifestVars=c("ob"),
                             mxPath("one", "ob", value=cl, free=FALSE),
                             mxPath("ob", arrows=2, value=noise, free=FALSE),
                             mxFitFunctionML(vector=TRUE))
  }

  m1 <-  
    mxModel("hmm", classes,
            mxData(data.frame(ob=trailN), "raw"),
            mxMatrix(nrow=3, ncol=1,
                     labels=paste0('i',1:3), name="initial"),
            mxMatrix(nrow=length(classes), ncol=length(classes),
                     labels=paste0('t', 1:(length(classes) * length(classes))),
                     name="transition"),
            mxExpectationHiddenMarkov(
              components=sapply(classes, function(m) m$name),
              initial="initial",
              transition="transition", scale="softmax"),
            mxFitFunctionML())

  m1$transition$free[1:(length(classes)-1), 1:length(classes)] <- TRUE

m1 <- mxRun(m1)

summary(m1)

print(m1$expectation$output)
}
