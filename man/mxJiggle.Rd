\name{mxJiggle}
\alias{mxJiggle}
\alias{imxJiggle}

\title{
Jiggle parameter values.
}

\description{
Jiggle free parameter values, subject to box constraints.  \code{imxJiggle()} is called internally by \code{\link{mxTryHard}()} (q.v.).  \code{mxJiggle()} provides a more user-friendly wrapper to \code{imxJiggle()}, and can alternately emulate the 'JIGGLE' behavior of classic Mx.
}

\usage{
mxJiggle(model, classic=FALSE, dsn=c("runif","rnorm","rcauchy"), loc=1, scale=0.25)
imxJiggle(params, lbounds, ubounds, dsn, loc, scale)
}

\arguments{
  \item{model}{An object of class MxModel.}
  \item{classic}{Logical; should \code{mxJiggle()} emulate the classic-Mx behavior elicited by keyword JIGGLE?  Defaults to \code{FALSE}.  See below, under "Details," for additional information.}
  \item{dsn}{Character string naming which random-number distribution--either \link[stats:runif]{uniform (rectangular)}, \link[stats:rnorm]{normal (Gaussian)}, or \link[stats:rcauchy]{Cauchy}--to be used to perturb free-parameter values.  Defaults to the uniform distribution (for \code{mxJiggle()}).}
   \item{loc, scale}{Numeric. The location and scale parameters of the distribution from which random values are drawn to perturb free-parameter values, defaulting respectively to 1 and 0.25 (for \code{mxJiggle()}).}
  \item{params}{Numeric vector of current free parameter values.}
  \item{lbounds}{Numeric vector of lower bounds on parameters.}
	\item{ubounds}{Numeric vector of upper bounds on parameters.}
}

\details{
If \code{mxJiggle()} argument \code{classic=FALSE} (the default), \code{mxJiggle()} calls \code{imxJiggle()}.  In that case, \code{mxJiggle()} passes \code{imxJiggle()} its own values for arguments \code{dsn}, \code{loc}, and \code{scale}, and extracts values for arguments \code{params}, \code{lbounds}, and \code{ubounds} from \code{model}.  Then, \code{model}'s free-parameter values are randomly perturbed before being re-assigned to it.  The distributional family from which the perturbations are randomly generated is dictated by argument \code{dsn}.  The distribution is parameterized by arguments \code{loc} and \code{scale}, respectively the location and scale parameters.  The location parameter is the distribution's median.  For the uniform distribution, \code{scale} is the absolute difference between its median and extrema (i.e., half the width of the rectangle); for the normal distribution, \code{scale} is its standard deviation; and for the Cauchy, \code{scale} is one-half its interquartile range.  Free-parameter values are first multiplied by random draws from a distribution with the provided \code{loc} and \code{scale}, then added to random draws from a distribution with the same \code{scale} but with a median of zero.

If \code{mxJiggle()} argument \code{classic=TRUE}, then each free-parameter value \eqn{x_i} is replaced with \eqn{x_i + 0.1(x_i + 0.5)}; this is the same behavior elicited in classic Mx by keyword JIGGLE.
}
\value{
\code{imxJiggle()} returns a numeric vector of randomly perturbed free-parameter values.  \code{mxJiggle()} returns \code{model}, with its free parameter values altered according to the other function arguments.
}

\seealso{
\code{\link{mxTryHard}()}
}
\examples{
data(demoOneFactor)
manifests <- names(demoOneFactor)
latents <- c("G")
factorModel <- mxModel(
	"One Factor",
	type="RAM",
	manifestVars = manifests,
	latentVars = latents,
	mxPath(from=latents, to=manifests,values=0.8),
	mxPath(from=manifests, arrows=2,values=1),
	mxPath(from=latents, arrows=2,
				 free=FALSE, values=1.0),
	mxData(cov(demoOneFactor), type="cov",
				 numObs=500)
)

iniPars <- coef(factorModel)
print(iniPars)

pars2 <- imxJiggle(params=iniPars,lbounds=NA,ubounds=NA,dsn="runif",loc=1,scale=0.05)
print(pars2)

mod2 <- mxJiggle(model=factorModel,scale=0.05)
coef(mod2)

mod3 <- mxJiggle(model=factorModel,classic=TRUE)
coef(mod3)
}
