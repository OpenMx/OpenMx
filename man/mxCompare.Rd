\name{mxCompare}
\alias{mxCompare}

\title{Assign Model Parameters}

\description{
    Compare the fit of one or more models to a base model.
    The output is a table with one row per model comparison, and the following columns:
\describe{
\item{base}{Name of the base model}
\item{comparison}{Name of the comparison model}
\item{ep}{Estimated parameters of the comparison model}
\item{minus2LL}{Minus 2*log-likelihood of the comparison model}
\item{df}{Degrees in freedom of the comparison model}
\item{AIC}{Akaike's Information Criterion for the comparison model}
\item{diffLL}{Change in minus 2*log-likelihood}
\item{diffdf}{Change in degrees of freedom}
\item{p}{Significance level of the change in fitness function}}
}

\usage{
mxCompare(base, comparison, digits = 3, all = FALSE)
}

\arguments{
   \item{base}{an MxModel object or list of MxModel objects.}
   \item{comparison}{a MxModel object or list of MxModel objects.}
   \item{digits}{a numeric setting table decimal precision.}
   \item{all}{A boolean value on whether to compare all bases with all comparisons.}
}

\seealso{
\code{\link{mxModel}} 
}

\examples{

data(demoOneFactor)
manifests <- names(demoOneFactor)
latents <- c("G1")
model1 <- mxModel("One Factor", type="RAM",
      manifestVars = manifests,
      latentVars = latents,
      mxPath(from=latents, to=manifests),
      mxPath(from=manifests, arrows=2),
      mxPath(from=latents, arrows=2, free=FALSE, values=1.0),
      mxData(cov(demoOneFactor), type="cov",numObs=500)
)

fit1 <- mxRun(model1)

latents <- c("G1","G2")
model2 <- mxModel("Two Factor", type="RAM",
      manifestVars = manifests,
      latentVars = latents,
      mxPath(from=latents[1], to=manifests[1:3]),
      mxPath(from=latents[2], to=manifests[4:5]),
      mxPath(from=manifests, arrows=2),
      mxPath(from=latents, arrows=2, free=FALSE, values=1.0),
      mxData(cov(demoOneFactor), type="cov",numObs=500)
)
fit2 <- mxRun(model2)

mxCompare(fit1, c(fit1, fit2), digits = 3)
}
