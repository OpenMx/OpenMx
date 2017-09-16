\name{mxExpectationHiddenMarkov}
\alias{MxExpectationHiddenMarkov-class}
\alias{mxExpectationHiddenMarkov}
\alias{print,MxExpectationHiddenMarkov-method}
\alias{show,MxExpectationHiddenMarkov-method}
\title{Hidden Markov expectation}
\usage{
mxExpectationHiddenMarkov(components, initial="initial", transition=NULL,
		      ..., verbose=0L, scale=c('softmax','sum'))
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
column sum.}
}
\description{
  Used in conjuction with \link{mxFitFunctionML}, this expectation
  can express a mixture model (with the transition matrix omitted)
  or a Hidden Markov model.
}
\details{
  Parameters are estimated in the given scale. To obtain the initial
  column vector and left stochastic transition matrix in probability
  units then examine the expectation's \code{output} slot.

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
