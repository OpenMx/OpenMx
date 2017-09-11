library(OpenMx)

data(demoOneFactor)
manifests <- names(demoOneFactor)
latents <- c("G")

template <- mxModel(
    "template", type="RAM",
    manifestVars = manifests,
    latentVars = latents,
    mxPath(from=latents, to=manifests, values=rnorm(length(manifests))),
    mxPath(from=manifests, arrows=2, values=rlnorm(length(manifests))),
    mxPath(from=latents, arrows=2, free=FALSE, values=1.0),
    mxPath(from = 'one', to = manifests, values=rnorm(length(manifests))))

factorRaw <- mxModel(template, name="OneFactorRaw",
		     mxData(demoOneFactor, type="raw"))

factorCov <- mxModel(template, name="OneFactorCov", 
		     mxData(observed=cov(demoOneFactor), means=colMeans(demoOneFactor),
			    type="cov", numObs=nrow(demoOneFactor)))

plan <- mxComputeSequence(list(
    mxComputeOnce('fitfunction', 'fit'),
    mxComputeNumericDeriv(checkGradient=FALSE, hessian=FALSE, iterations=2),
    mxComputeReportDeriv(),
    mxComputeReportExpectation()
))

factorRaw <- mxRun(mxModel(factorRaw, plan))
factorCov <- mxRun(mxModel(factorCov, plan))

omxCheckCloseEnough(factorRaw$output$fit, factorCov$output$fit + prod(dim(demoOneFactor))*log(2*pi), 1e-10)
omxCheckCloseEnough(factorRaw$output$gradient, factorCov$output$gradient, 1e-5)
