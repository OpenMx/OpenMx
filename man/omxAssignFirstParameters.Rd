\name{omxAssignFirstParameters}
\alias{omxAssignFirstParameters}

\title{Assign First Available Values to Model Parameters}

\description{
    Assign starting values to the free parameters of a model. Select one of the current values for each
    free parameter and use that value.
}

\usage{
omxAssignFirstParameters(model, indep = FALSE)
}

\arguments{
   \item{model}{a MxModel object.}
   \item{indep}{assign parameters to independent submodels.}
}

\seealso{
\code{\link{omxGetParameters}}, \code{\link{omxSetParameters}} 
}

\examples{

A <- mxMatrix('Full', 3, 3, labels = c('a','b', NA), free = TRUE, name = 'A')
model <- mxModel(A, name = 'model')
model <- omxAssignFirstParameters(model)

}
