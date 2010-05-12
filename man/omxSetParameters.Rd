\name{omxSetParameters}
\alias{omxSetParameters}

\title{Assign Model Parameters}

\description{
    Modify the attributes of parameters in a model.  This function cannot modify parameters that have NA labels.
}

\usage{
omxSetParameters(model, labels, free = NULL, values = NULL,
	newlabels = NULL, lbound = NULL, ubound = NULL, indep = FALSE)
}

\arguments{
   \item{model}{a MxModel object.}
   \item{labels}{a character vector of target parameter names.}
   \item{free}{a boolean vector of parameter free/fixed designations.} 
   \item{values}{a numeric vector of parameter values.}
   \item{newlabels}{a character vector of new parameter names.} 
   \item{lbound}{a numeric vector of lower bound values.} 
   \item{ubound}{a numeric vector of upper bound values.} 
   \item{indep}{boolean. set parameters in independent submodels.}
}

\seealso{
\code{\link{omxGetParameters}}, \code{\link{omxAssignFirstParameters}} 
}

\examples{

A <- mxMatrix('Full', 3, 3, labels = c('a','b', NA), free = TRUE, name = 'A')
model <- mxModel(A, name = 'model')
model <- omxSetParameters(model, c('a', 'b'), values = c(1, 2))
model <- omxSetParameters(model, c('a', 'b'), newlabels = c('b', 'a'))

}
