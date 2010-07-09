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

A <- mxMatrix('Full', 2, 2, labels = c("A11", "A12", "A21", NA), values= 1:4, free = c(TRUE,TRUE,FALSE,TRUE), byrow=TRUE, name = 'A')
model <- mxModel(A, name = 'model')
omxGetParameters(model)

# A11  A12 <NA> 
#   1    2    4 

A@labels
#      [,1]  [,2] 
# [1,] "A11" "A12"
# [2,] "A21" NA   

A@free
#       [,1] [,2]
# [1,]  TRUE TRUE
# [2,] FALSE TRUE

A@labels
#      [,1]  [,2] 
# [1,] "A11" "A12"
# [2,] "A21" NA   

}
