\name{omxAssignFirstParameters}
\alias{omxAssignFirstParameters}

\title{Assign First Available Values to Model Parameters}

\description{
    Sometimes you may have a free parameter with two different starting values in your model, and it may be sufficient to arbitrarily select one of those starting values for optimization.  OpenMx will not run a model until all instances of a free parameter have the same starting value.  This function will assign starting values to the free parameters of a model. It will select an arbitrary current value (the "first" value it finds, where "first" is not defined) for each free parameter and use that value.
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

A <- mxMatrix('Full', 3, 3, values = c(1:9), 
   	labels = c('a','b', NA), free = TRUE, name = 'A')
model <- mxModel(A, name = 'model')
model <- omxAssignFirstParameters(model)

}
