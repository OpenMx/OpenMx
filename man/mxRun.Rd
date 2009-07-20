\name{mxRun}
\alias{mxRun}

\title{Send a Model to the Optimizer}

\description{
   This function begins optimization on the top-level model.
}

\usage{
mxRun(model)
}

\arguments{
    \item{model}{An \link{MxModel} object to be optimized.}
}

\details{
The mxRun function is used to optimize free parameters in \link{MxModel} objects based on an objective function. MxModel objects included in the mxRun function must include an appropriate objective function.

Free parameters are estimated or updated based on the objective function. These estimated values, along with estimation information and model fit, can be found in the 'output' slot of MxModel objects after mxRun has been used.

If a model is dependent on or shares parameters with another model, both models must be included as arguments in another MxModel object. This top-level MxModel object must include objective functions in both submodels, as well as an additional objective function describing how the results of the first two should be combined.
}

\value{
   Returns an MxModel object with free parameters updated to their final values.
   The return value contains an "output" slot with the results of 
   optimization.
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
