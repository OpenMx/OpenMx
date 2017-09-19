library(OpenMx)

omxCheckError(vcov(mxModel("empty")),
	      'This model has not been run yet. Tip: Use
  model = mxRun(model)
to estimate a model.')

data(demoOneFactor)
manifests <- names(demoOneFactor)
latents <- c("G")

oneFactor <- mxModel(
    "oneFactor", type="RAM",
    manifestVars = manifests,
    latentVars = latents,
    mxPath(from=latents, to=manifests, values=rnorm(length(manifests))),
    mxPath(from=manifests, arrows=2, values=rlnorm(length(manifests))),
    mxPath(from=latents, arrows=2, free=FALSE, values=1.0),
    mxPath(from = 'one', to = manifests, values=rnorm(length(manifests))),
    mxData(demoOneFactor, type="raw"))

oneFactor <- mxOption(oneFactor, "Calculate Hessian", "No")

oneFactor <- mxRun(oneFactor)

if (mxOption(NULL, "Default optimizer") != 'SLSQP') {
  vc1 <- omxCheckWarning(vcov(oneFactor),
                         "The 'Calculate Hessian' option is disabled. This may result in a poor accuracy vcov matrix.
Turn on with mxOption(model, 'Calculate Hessian', 'Yes')")
}

oneFactor <- mxOption(oneFactor, "Calculate Hessian", "Yes")

oneFactor <- mxRun(oneFactor)

omxCheckTrue(all(abs(sqrt(diag(vcov(oneFactor))) - oneFactor$output$standardErrors) < 1e-6))
