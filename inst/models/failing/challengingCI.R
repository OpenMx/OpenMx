library(OpenMx)

if (0) {
  # move parameters to MLE
  load("memprobmodel2.RData")
  memprobmodel2 <- mxRun(memprobmodel2, intervals=FALSE)
  save(memprobmodel2, file="memprobmodel2.RData")
}

got <- try(load("models/failing/memprobmodel2.RData"))
if (is(got, "try-error")) {
  load("memprobmodel2.RData")
}

fit <- mxRun(mxModel(memprobmodel2, mxComputeOnce('fitfunction','fit')))
target <- fit$output$fit + 3.841459

for (engine in c("NPSOL", "CSOLNP")) {
  ciModel <- mxModel("ciModel", memprobmodel2,
                     mxAlgebra((target - ctsem.objective)^2 - ctsem.DRIFT[2,1], name="ci"),
                     mxFitFunctionAlgebra("ci"),
                     mxComputeGradientDescent(engine=engine))
  fit <- mxRun(ciModel)
  omxCheckCloseEnough(fit$ctsem$DRIFT$result[2,1], .3880241, .001)
}
