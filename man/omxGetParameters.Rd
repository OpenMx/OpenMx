\name{omxGetParameters}
\alias{omxGetParameters}

\title{Fetch Model Parameters}

\description{
    Return a vector of the free parameters in the model.
}

\usage{
omxGetParameters(model, indep = FALSE)
}

\arguments{
   \item{model}{a MxModel object}
   \item{indep}{fetch parameters from independent submodels.}
}

\seealso{
\code{\link{omxSetParameters}}, \code{\link{omxAssignFirstParameters}} 
}

\examples{

A <- mxMatrix('Full', 3, 3, labels = c('a','b', NA), free = TRUE, name = 'A')
model <- mxModel(A, name = 'model')
parameters <- omxGetParameters(model)

}
