library(OpenMx)

set.seed(1)
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

factorCov <- mxModel(template, name="OneFactorCov", 
		     mxData(observed=cov(demoOneFactor), means=colMeans(demoOneFactor),
			    type="cov", numObs=nrow(demoOneFactor)))

orig <- coef(factorCov)

fit <- rep(NA, length(orig))
for (vx in 1:length(orig)) {
  mask <- rep(TRUE, length(orig))
  mask[vx] <- FALSE
  m1 <- omxSetParameters(factorCov, labels=names(orig), free=mask)
  m1 <- mxRun(mxModel(m1, mxComputeGradientDescent()), silent = TRUE)
	
  plan <- mxComputeGradientDescent()
  plan$.excludeVars <- names(orig)[vx]
  
  m2 <- mxRun(mxModel(factorCov, plan), silent = TRUE)
  omxCheckCloseEnough(m1$output$fit, m2$output$fit, 1e-4)
  fit[vx] <- m1$output$fit
}

omxCheckCloseEnough(log(sd(fit)), 6.3, .5)
