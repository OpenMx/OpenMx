\name{mxCompare}
\alias{mxCompare}

\title{Assign Model Parameters}

\description{
    Compare the fit of one or more models to a base model.
    The output is a table with one row per model comparison, and the following columns: 
         base  - Name of the base model
         name  - Name of the model model being compared on this line
         ep    - Estimated parameters of the comparison model
         -2LL  - minus 2*log-likelihood of the current model
         df    - Change in degrees of freedom
        AIC    - Akaike's Information Criterion for the comparison
       diffLL  - Change in minus 2* log-likelihood
       diffdf  - Change in degrees of freedom
       p       - Significance level of the change in fit
}

\usage{
omxCompare(base,comparison,digits=3)
}

\arguments{
   \item{base}{an MxModel object.}
   \item{comparison}{a vector of one or more nested MxModel objects}
   \item{digits}{a numeric setting table decimal precision.} 
}

\seealso{
\code{\link{mxModel}} 
}

\examples{

require(OpenMx)
data(demoOneFactor)
manifests <- names(demoOneFactor)
latents <- c("G1")
model1 <- mxModel("One Factor", type="RAM",
      manifestVars = manifests,
      latentVars = latents,
      mxPath(from=latents, to=manifests),
      mxPath(from=manifests, arrows=2),
      mxPath(from=latents, arrows=2, free=F, values=1.0),
      mxData(cov(demoOneFactor), type="cov",numObs=500)
)
fit1 = mxRun(model1)

latents <- c("G1","G2")
model2 <- mxModel("Two Factor", type="RAM",
      manifestVars = manifests,
      latentVars = latents,
      mxPath(from=latents[1], to=manifests[1:3]),
      mxPath(from=latents[2], to=manifests[4:5]),
      mxPath(from=manifests, arrows=2),
      mxPath(from=latents, arrows=2, free=F, values=1.0),
      mxData(cov(demoOneFactor), type="cov",numObs=500)
)
fit2 = mxRun(model2)

mxCompare(fit1, c(fit1,fit2), digits=3)
#         base comparison ep minus2LL df   AIC diffLL diffdf p
# 1 One Factor One Factor 10    -3650  5 -2.62      0      0 0
# 2 One Factor Two Factor 10    -2450  5  1200   1200      0 0
}
