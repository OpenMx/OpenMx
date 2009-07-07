\name{mxOption}
\alias{mxOption}

\title{Set or Clear an Optimizer Option}

\description{
The function sets or clears an option that is specific to the optimizer in the back-end.
}

\usage{
mxOption(model, key, value)
}

\arguments{
\item{model}{An \link{MxModel} object}
\item{key}{The name of the option.}
\item{value}{The value of the option.}
}

\details{
Sets an option that is specific to the particular optimizer used in the back-end. The name of the option is the \sQuote{key} argument. Use value = NULL to remove an existing option. Before the model is submitted to the back-end, all keys and values are converted into strings using the \link{as.character} function. Use \code{getOption('mxOptimizerOptions')} to see the default optimizer options.
}
	
\value{
Returns the model with the optimizer option either set or cleared.
}

\references{
More information about Mx may be found at http://www.vcu.edu/mx/documentation.html
}


\examples{
model <- mxModel()

model <- mxOption(model, "Function Precision", 1e-5)
model <- mxOption(model, "Hessian", "No")
model <- mxOption(model, "Function Precision", NULL)

}
