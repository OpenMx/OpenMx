\name{omxSetParameters}
\alias{omxSetParameters}

\title{Assign Model Parameters}

\description{
    Assign starting values to the free parameters of a model.  This function cannot assign values to free parameters with NA labels.
}

\usage{
omxSetParameters(model, names, values, indep = FALSE)
}

\arguments{
   \item{model}{a MxModel object.}
   \item{names}{a vector of free parameter names.}
   \item{values}{a vector of free parameter values.}   
   \item{indep}{set parameters in independent submodels.}
}

\seealso{
\code{\link{omxGetParameters}}, \code{\link{omxAssignFirstParameters}} 
}

\examples{

A <- mxMatrix('Full', 3, 3, labels = c('a','b', NA), free = TRUE, name = 'A')
model <- mxModel(A, name = 'model')
model <- omxSetParameters(model, c('a', 'b'), c(1, 2))

}
