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
   \item{free}{fetch either free parameters (TRUE), or fixed parameters or both types. Default value is TRUE.}
}

\details{
  The argument \sQuote{free} dictates whether to return only free parameters or only fixed parameters or both free and fixed parameters. The function can return unlabelled free parameters (parameters with a label of NA). These anonymous free parameters will be identified as \sQuote{modelname.matrixname[row,col]}.  It will not return fixed parameters that have a label of NA. No distinction is made between ordinary labels, definition variables, and square bracket constraints.
}

\seealso{
\code{\link{omxSetParameters}}, \code{\link{omxLocateParameters}}, \code{\link{omxAssignFirstParameters}} 
}

\examples{

A <- mxMatrix('Full', 2, 2, labels = c("A11", "A12", "A21", NA), values= 1:4, free = c(TRUE,TRUE,FALSE,TRUE), byrow=TRUE, name = 'A')
model <- mxModel(A, name = 'model')

# Request all free parameters in model
omxGetParameters(model)

# A11  A12 model.A[2,2] 
#   1    2    4 

# Request fixed parameters from model
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

# Example using un-labelled parameters
# Load the OpenMx Package if not already loaded
require(OpenMx) 
# Read in some demo data
data(demoOneFactor)
# Grab the names for manifestVars 
manifestVars <- names(demoOneFactor)
nVar = length(manifestVars) # 5 variables
factorModel <- mxModel("One Factor",
    mxMatrix(name="A", type="Full", nrow=nVar, ncol=1, values=0.2, free=TRUE, labels=letters[1:nVar]),
    mxMatrix(name="L", type="Symm", nrow=1, ncol=1, values=1, free=FALSE),
    # the "U" matrix has nVar (5) anonymous free parameters
    mxMatrix(name="U", type="Diag", nrow=nVar, ncol=nVar, values=1, free=TRUE),
    mxAlgebra(expression=A \%&\% L + U, name="R"),
    mxMLObjective(covariance="R", dimnames=manifestVars),
    mxData(observed=cov(demoOneFactor), type="cov", numObs=500)
)

# Get all free parameters
params         <- omxGetParameters(factorModel)
# Set new values for these params, saving them in a new model
newFactorModel <- omxSetParameters(factorModel, names(params), values = 1:10)
# Read out the values from the new model
newParams      <- omxGetParameters(newFactorModel)

}
