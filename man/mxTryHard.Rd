\name{mxTryHard}
\alias{mxTryHard}

\title{Make multiple attempts to run a model}

\description{
Makes multiple attempts to fit an MxModel object with \code{\link{mxRun}()} until the optimizer yields an acceptable solution or the maximum number of attempts is reached.  Each attempt uses the parameter estimates of the previous attempt as start values, but they are each multiplied by random draws from a uniform distribution, and optimization parameters may be altered.  From among its attempts, the function returns the fitted, post-\code{\link{mxRun}()} model with the smallest fit-function value, and can print to the console the start values it used for that model.
}

\usage{
mxTryHard(model, extraTries = 10, greenOK = FALSE, loc = 1, 
  scale = 0.25, initialGradientStepSize = .00001, initialGradientIterations = 1,
    initialTolerance=1e-12, checkHess = TRUE, fit2beat = Inf, paste = TRUE,
  iterationSummary=FALSE, bestInitsOutput=TRUE, showInits=FALSE, verbose = 0, intervals = FALSE)
}

\arguments{
  \item{model}{The model to be run; object of class \code{\link{MxModel}}.}
  \item{extraTries}{The number of attempts to run the model \emph{in addition to the first}.  In effect, is the maximum number of attempts \code{mxTryHard()} will make, since the function will stop once an acceptable solution is reached.  Defaults to 10, in which case a maximum of 11 total attempts will be made.}
  \item{greenOK}{Logical; is a solution with Mx status GREEN (npsolstatus=1) acceptable?  Defaults to \code{FALSE}.}
  \item{loc, scale}{The location and scale parameters of the uniform (rectangular) distribution from which random values are drawn to disturb start values between attempts.  The location parameter is the distribution's median, and the scale parameter is the half-width of the rectangle (that is, the absolute difference between the median and the extrema).  Defaults to a uniform distribution on the interval (0.75, 1.25).}
  \item{initialGradientStepSize, initialGradientIterations, initialTolerance}{Optimization parameters passed to \code{\link{mxComputeGradientDescent}}.}
  \item{checkHess}{Logical; is a positive-definite Hessian a requirement for an acceptable solution?  Defaults to \code{TRUE}.}
  \item{fit2beat}{An upper limit to the objective-function value that an acceptable solution may have.  Useful if a nested submodel of \code{model} has already been fitted, since \code{model}, with its additional free parameters, should not yield a fit-function value any greater than that of the submodel.}
  \item{paste}{Logical.  If \code{TRUE} (default), start values for the returned fitted model are printed to console as a comma-separated string.  This is useful if the user wants to copy-paste these values into an \R script, say, in an \code{\link{omxSetParameters}()} statement. If \code{FALSE}, the vector of start values is printed as-is.  Note that this vector, from \code{\link{omxGetParameters}()}, has names corresponding to the free parameters; these names are not displayed when \code{paste=TRUE}.}
  \item{iterationSummary}{Logical. If \code{TRUE}, displays parameter estimates and fit values for every fit attempt. Defaults to \code{FALSE}.}
  \item{bestInitsOutput}{Logical. If \code{TRUE}, outputs starting values that resulted in best fit, according to format specified by \code{paste} argument. Defaults to \code{TRUE}.}
  \item{showInits}{Logical. If \code{TRUE}, displays starting values for every fit attempt. Defaults to \code{FALSE}.}
  \item{verbose}{Passed to \code{\link{mxComputeGradientDescent}} to specify level of output to console during optimization.}
  \item{intervals}{Logical. If TRUE, OpenMx will estimate any specified confidence intervals.}
}


\value{
Usually, \code{mxTryHard()} returns a post-\code{\link{mxRun}()} \code{\link{MxModel}} object.  Specifically, this will be the fitted model having the smallest fit-function value found by \code{mxTryHard()} during its attempts.  The start values used to obtain this fitted model are printed to console.  

If every attempt at running \code{model} fails, \code{mxTryHard()} returns an object of class 'try-error', and the start values from the last attempt are printed to console.

\code{mxTryHard()} throws a warning if the returned \code{\link{MxModel}} object has a nonzero npsolstatus.
}

\seealso{
\code{\link{mxRun}()}
}

\examples{ %Could use a better example...
library(OpenMx)

data(demoOneFactor)  # load the demoOneFactor dataframe

manifests <- names(demoOneFactor) # set the manifest to the 5 demo variables
latents   <- c("G")  # define 1 latent variable
model <- mxModel(model="One Factor", type="RAM",
    manifestVars = manifests,
    latentVars   = latents,
    mxPath(from=latents, to=manifests, labels=paste("b", 1:5, sep="")),
    mxPath(from=manifests, arrows=2, labels=paste("u", 1:5, sep="")),
    mxPath(from=latents  , arrows=2, free=FALSE, values=1.0),
    mxData(cov(demoOneFactor), type="cov", numObs=500)
)
model <- mxTryHard(model) # Run the model, returning the result into model
summary(model) # Show summary of the fitted model

}
