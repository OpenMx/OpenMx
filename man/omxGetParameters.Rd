\name{omxGetParameters}
\alias{omxGetParameters}

\title{Fetch Model Parameters}

\description{
    Return a vector of the free parameters in the model.
}

\usage{
omxGetParameters(model, indep = FALSE, free = c(TRUE, FALSE, NA))
}

\arguments{
   \item{model}{a MxModel object}
   \item{indep}{fetch parameters from independent submodels.}
   \item{free}{fetch either free parameters or fixed parameters or both types.}
}

\details{
  The argument \sQuote{free} dictates whether to return only free parameters or only fixed parameters or both free and fixed parameters. The function will return free parameters that have a label of NA.  But it will never return fixed parameters that have a label of NA. No distinction is made between ordinary labels, and definition variables, and square bracket constraints in labels.
}

\seealso{
\code{\link{omxSetParameters}}, \code{\link{omxAssignFirstParameters}} 
}

\examples{

A <- mxMatrix('Full', 2, 2, labels = c("A11", "A12", "A21", NA), values= 1:4, free = c(TRUE,TRUE,FALSE,TRUE), byrow=TRUE, name = 'A')
model <- mxModel(A, name = 'model')

# Request all free paraemters in model
omxGetParameters(model)

# A11  A12 <NA> 
#   1    2    4 

# Request fixed paraemters from model
omxGetParameters(model, free=FALSE)
# A21 
#   3

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
