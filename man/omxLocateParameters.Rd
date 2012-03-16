\name{omxLocateParameters}
\alias{omxLocateParameters}

\title{Summarize Model Parameters}

\description{
    Return a data.frame object summarizing the free parameters in the model.
}

\usage{
omxLocateParameters(model, labels = NULL, indep = FALSE)
}

\arguments{
   \item{model}{a MxModel object}
   \item{labels}{optionally specify which free parameters to retrieve.}
   \item{indep}{fetch parameters from independent submodels.}
}

\details{
   Invoking the function with the default value for the \sQuote{labels}
   argument retrieves all the free parameters.  The \sQuote{labels}
   argument can be used to select a subset of the free parameters.
   Note that \sQuote{NA} is a valid possible input to the \sQuote{labels} argument.
}

\seealso{
\code{\link{omxGetParameters}}, \code{\link{omxSetParameters}}, \code{\link{omxAssignFirstParameters}} 
}

\examples{

A <- mxMatrix('Full', 2, 2, labels = c("A11", "A12", NA, NA), values= 1:4,
   free = TRUE, byrow = TRUE, name = 'A')

model <- mxModel(A, name = 'model')

# Request all free parameters in model
omxLocateParameters(model)

# Request free parameters "A11" and all NAs
omxLocateParameters(model, c("A11", NA))


}
