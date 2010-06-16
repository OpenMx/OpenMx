require(OpenMx)

data(demoOneFactor)

manifests <- names(demoOneFactor)
latents <- c("factor")

factorModel <- mxModel("One Factor", type="RAM",
    manifestVars = manifests,
    latentVars = latents,
    mxPath(from=latents, to=manifests),
    mxPath(from=manifests, arrows=2, lbound=.0001),
    mxPath(from=latents, arrows=2, free=F, values=1.0),
    mxData(cov(demoOneFactor), type="cov", numObs=500),
    mxAlgebra(sqrt(diag2vec(S)),name="P"),
    mxCI(c("A", "S", "P"))
)
factorFitCI <- mxRun(factorModel, intervals=TRUE)
summary(factorFitCI)