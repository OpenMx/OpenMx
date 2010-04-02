\name{mxRun}
\alias{mxRun}

\title{Send a Model to the Optimizer}

\description{
   This function begins optimization on the top-level model.
}

\usage{
mxRun(model, ..., intervals = FALSE, silent = FALSE, unsafe = FALSE)
}

\arguments{
    \item{model}{A \link{MxModel} object to be optimized.}    
    \item{...}{Not used.  Forces remaining arguments to be specified by name.}
    \item{intervals}{A boolean indicating whether to compute the specified confidence intervals.}
    \item{silent}{A boolean indicating whether to print status to terminal.}
    \item{unsafe}{A boolean indicating whether to ignore errors.}
}

\details{
The mxRun function is used to optimize free parameters in \link{MxModel} objects based on an objective function. MxModel objects included in the mxRun function must include an appropriate objective function.

If the \sQuote{silent} flag is TRUE, then model execution will not print any status messages to the terminal, and a warning is not issued if NPSOL returns a non-zero status code.

If the \sQuote{unsafe} flag is TRUE, then any error conditions will throw a warning instead of an error.  It is strongly recommended to use this feature only for debugging purposes.

Free parameters are estimated or updated based on the objective function. These estimated values, along with estimation information and model fit, can be found in the 'output' slot of MxModel objects after mxRun has been used.

If a model is dependent on or shares parameters with another model, both models must be included as arguments in another MxModel object. This top-level MxModel object must include objective functions in both submodels, as well as an additional objective function describing how the results of the first two should be combined.
}

\value{
   Returns an MxModel object with free parameters updated to their final values.
   The return value contains an "output" slot with the results of 
   optimization.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{
#Create a model that includes data, matrices A, S and F, and an objective function
\dontrun{
data <- mxData(mydata, type="cov", numObs = 100)
objective <- mxRAMObjective('A', 'S', 'F')
model <- mxModel("mymodel", A, S, F, data, objective)

#Use mxRun to optimize the free parameters in the matrices A and S
model <- mxRun(model)

#print the output
model@output
}
}
