\name{mxTryHard}
\alias{mxTryHard}
\alias{mxTryHardOrig}
\alias{mxTryHardctsem}
\alias{mxTryHardWideSearch}
\alias{mxTryHardOrdinal}
\alias{THard}

\title{Make multiple attempts to run a model}

\description{
Makes multiple attempts to fit an MxModel object with \code{\link{mxRun}()} until the optimizer yields an acceptable solution or the maximum number of attempts (set by \code{extraTries}) is reached.  Between attempts, start values are perturbed by random numbers (see details). Optimization-control parameters may also be altered.  From among its attempts, the function returns the fitted model with the best fit (smallest fit-function value). If \code{bestInitsOutput} is TRUE, then in addition the start values used for the best-fitting model will be printed to the console.

\emph{note}: If the model contains \code{\link{mxConstraint}}s, or if \code{\link{mxFitFunctionWLS}} is being used, the Hessian cannot be checked, and so \code{checkHess} will be coerced to FALSE.
}

\usage{
mxTryHard(model, extraTries = 10, greenOK = FALSE, loc = 1, scale = 0.25,
initialGradientStepSize = imxAutoOptionValue("Gradient step size"),
initialGradientIterations = imxAutoOptionValue('Gradient iterations'),
initialTolerance=as.numeric(mxOption(NULL,'Optimality tolerance')), 
checkHess = TRUE, fit2beat = Inf, paste = TRUE,	iterationSummary=FALSE, 
bestInitsOutput=TRUE, showInits=FALSE, verbose=0, intervals = FALSE,
finetuneGradient=TRUE, jitterDistrib=c("runif","rnorm","rcauchy"), exhaustive=FALSE,
maxMajorIter=3000, OKstatuscodes, wtgcsv=c("prev","best","initial"), silent=interactive())

mxTryHardOrig(model, finetuneGradient=FALSE, maxMajorIter=NA, 
wtgcsv=c("prev","best"), silent=FALSE, ...)

mxTryHardctsem(model, initialGradientStepSize = .00001, 
initialGradientIterations = 1,
initialTolerance=1e-12,	jitterDistrib="rnorm", ...)

mxTryHardWideSearch(model, finetuneGradient=FALSE, jitterDistrib="rcauchy", 
exhaustive=TRUE, wtgcsv="prev", ...)

mxTryHardOrdinal(model, greenOK = TRUE,	checkHess = FALSE, 
finetuneGradient=FALSE, exhaustive=TRUE,
OKstatuscodes=c(0,1,5,6), wtgcsv=c("prev","best"), ...)
}

\arguments{
  \item{model}{The MxModel to be run.}
  \item{extraTries}{The number of attempts to run the model \emph{in addition to the first}.  In effect, is the maximum number of attempts \code{mxTryHard()} will make, since the function will stop once an acceptable solution is reached.  Defaults to 10 (for \code{mxTryHard()}), in which case a maximum of 11 total attempts will be made.}
  \item{greenOK}{Logical; is a solution with Mx status GREEN (optimizer status code 1) acceptable?  Defaults to \code{FALSE} (for \code{mxTryHard()}).  Ignored if a value is provided for \code{OKstatuscodes}.}
  \item{loc, scale}{Numeric. The location and scale parameters of the distribution from which random values are drawn to perturb start values between attempts, defaulting respectively to 1 and 0.25.  See below, under "Details," for additional information.}
  \item{initialGradientStepSize, initialGradientIterations, initialTolerance}{Numeric.  Initial values of optimization-control parameters passed to \code{\link{mxComputeGradientDescent}()} if \code{model} is using the default compute plan.}
  \item{checkHess}{Logical; is a positive-definite Hessian a requirement for an acceptable solution?  Defaults to \code{TRUE} (for \code{mxTryHard()}).  If \code{TRUE}, the Hessian and standard errors are calculated with each fit attempt, irrespective of the value of relevant \link[=mxOption]{options}.  The exception is if \code{model} or any of its submodels contains \code{\link{MxConstraint}s}, in which case \code{checkHess} is coerced to \code{FALSE}. }
  \item{fit2beat}{Numeric upper limit to the fitfunction value that an acceptable solution may have.  Useful if a nested submodel of \code{model} has already been fitted, since \code{model}, with its additional free parameters, should not yield a fitfunction value any greater than that of the submodel.}
  \item{paste}{Logical.  If \code{TRUE} (default), start values for the returned fitted model are printed to console as a comma-separated string.  This is useful if the user wants to copy-paste these values into an \R script, say, in an \code{\link{omxSetParameters}()} statement. If \code{FALSE}, the vector of start values is printed as-is.  Note that this vector, from \code{\link{omxGetParameters}()}, has names corresponding to the free parameters; these names are not displayed when \code{paste=TRUE}.}
  \item{iterationSummary}{Logical. If \code{TRUE}, displays parameter estimates and fit values for every fit attempt, even if \code{silent=TRUE}. Defaults to \code{FALSE}.}
  \item{bestInitsOutput}{Logical. If \code{TRUE} and if \code{silent=FALSE}, \code{mxTryHard()} displays the starting values that resulted in the best fit, according to format specified by \code{paste} argument. Defaults to \code{TRUE}.}
  \item{showInits}{Logical. If \code{TRUE}, displays starting values for every fit attempt, even if \code{silent=TRUE}.  Defaults to \code{FALSE}.}
  \item{verbose}{If \code{model} is using the default compute plan, is passed to \code{\link{mxComputeGradientDescent}()} to specify level of output printed to console during optimization.}
  \item{intervals}{Logical. If TRUE, OpenMx will estimate any specified confidence intervals.}
  \item{finetuneGradient}{Logical.  If \code{TRUE} (default for \code{mxTryHard()}), then as repeated fit attempts appear to be improving, \code{mxTryHard()} will adjust optimization-control parameters \code{\link[OpenMx:mxComputeGradientDescent]{gradientStepSize}}, \code{\link[OpenMx:mxComputeGradientDescent]{gradientIterations}}, and \code{\link[OpenMx:mxComputeGradientDescent]{tolerance}}, as well as argument \code{scale}, to "fine-tune" its convergence toward an optimal solution.  \code{finetuneGradient=FALSE} is recommended for analyses involving thresholds.}
  \item{jitterDistrib}{Character string naming which random-number distribution--either \link[stats:runif]{uniform (rectangular)}, \link[stats:rnorm]{normal (Gaussian)}, or \link[stats:rcauchy]{Cauchy}--to be used to perturb start values.  Defaults to the uniform distribution (for \code{mxTryHard()}).  See below, under "Details," for additional information.}
  \item{exhaustive}{Logical.  If \code{FALSE} (default for \code{mxTryHard()}), \code{mxTryHard()} stops making additional attempts once it reaches an acceptable solution.  If \code{TRUE}, the function instead continues until it reaches its maximum number of attempts (as per \code{extraTries}), and returns the best solution it found.}
  \item{maxMajorIter}{Integer; passed to \code{\link{mxComputeGradientDescent}()}.  Defaults to 3000, which was the internally hard-coded value \code{mxTryHard()} used in a prior version of OpenMx.  Value of \code{NA} is permitted, in which case \code{mxTryHard()} will calculate a value via the on-load default formula for the "Major iterations" \link[=mxOption]{option}.}
  \item{OKstatuscodes}{Optional integer vector containing optimizer status codes that an acceptable solution is permitted to have.  \code{mxTryHard()} always considers a status code of 0 to be acceptable, this argument notwithstanding. By default, \code{mxTryHard()} will consider status code 0 acceptable, and, if \code{greenOK=TRUE}, status code 1 as well.  If a value is supplied for \code{OKstatuscodes} that conflicts with \code{greenOK}, \code{OKstatuscodes} controls.}
  \item{wtgcsv}{Character vector.  "Where to get current start values."  See below, under "Details," for additional information.}
  \item{silent}{Logical; for \code{mxTryHard()}, defaults to \code{TRUE} if running interactively, and to \code{FALSE} otherwise.  If \code{TRUE}, persistent \link[base:message]{message-printing} during execution of \code{mxTryHard()} is suppressed, and non-persistent printing is used instead. The two exceptions are the persistent printing requested by \code{TRUE} values of \code{iterationSummary} and \code{showInits}.}
  \item{...}{Additional arguments to be passed to \code{mxTryHard()}.}
}


\details{
\code{mxTryHardOrig()}, \code{mxTryHardctsem()}, \code{mxTryHardWideSearch()}, and \code{mxTryHardOrdinal()} are wrapper functions to the main workhorse function \code{mxTryHard()}.  Each wrapper function has default values for certain arguments that are tailored toward a specific purpose.  \code{mxTryHardOrig()} imitates the functionality of the earliest implementations of \code{mxTryHard()} in OpenMx's history; its chief purpose is to find good start values that lead to an acceptable solution.  \code{mxTryHardctsem()} uses \code{mxTryHard()} to "zero in" on an acceptable solution with models that can be difficult to optimize, such as continuous-time state-space models.  \code{mxTryHardWideSearch()} uses \code{mxTryHard()} to search a wide region of the parameter space, in hope of avoiding local fitfunction minima.  \code{mxTryHardOrdinal()} attempts to use \code{mxTryHard()} as well as it can be used with models involving ordinal data.

Argument \code{wtgcsv} dictates where \code{mxTryHard()} is permitted to find free-parameter values, at the start of each fit attempt after the first, before randomly perturbing them to create the current fit attempt's start values.  If \code{"prev"} is included, then \code{mxTryHard()} is permitted to use the parameter estimates of the most recent non-error fit attempt.  If \code{"best"} is included, then \code{mxTryHard()} is permitted to use the parameter estimates at the best solution so far.  If \code{"initial"} is included, then \code{mxTryHard()} is permitted to use the initial start values in \code{model}, as provided by the user.  The default is to permit all three, in which case \code{mxTryHard()} is written to use the best solution's values if available, and otherwise to use the most recent solution's values,  but to periodically revert to the initial values if recent fit attempts have not improved on the best solution.

Once the start values are located for the current fit attempt, they are randomly perturbed before being assigned to the MxModel.  The distributional family from which the perturbations are randomly generated is dictated by argument \code{jitterDistrib}.  The distribution is parameterized by arguments \code{loc} and \code{scale}, respectively the location and scale parameters.  The location parameter is the distribution's median.  For the uniform distribution, \code{scale} is the absolute difference between its median and extrema (i.e., half the width of the rectangle); for the normal distribution, \code{scale} is its standard deviation; and for the Cauchy, \code{scale} is one-half its interquartile range.  Start values are first multiplied by random draws from a distribution with the provided \code{loc} and \code{scale}, then added to random draws from a distribution with the same \code{scale} but with a median of zero.
}


\value{
Usually, \code{mxTryHard()} returns a post-\code{\link{mxRun}()} \code{\link{MxModel}} object.  Specifically, this will be the fitted model having the smallest fit-function value found by \code{mxTryHard()} during its attempts.  The start values used to obtain this fitted model are printed to console if \code{bestInitsOutput=TRUE}.  

If every attempt at running \code{model} fails, \code{mxTryHard()} returns an object of class 'try-error'.

\code{mxTryHard()} throws a warning if the returned \code{\link{MxModel}} object has a nonzero status code (unless nonzero status codes are considered acceptable per argument \code{greenOK} or \code{OKstatuscodes}).
}

\seealso{
\code{\link{mxRun}()}, \code{\link{mxComputeTryHard}}
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
