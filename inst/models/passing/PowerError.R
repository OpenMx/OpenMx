library(OpenMx)
library(testthat)
data(demoOneFactor)

latents  = c("G")
manifests = names(demoOneFactor)

m1 <- mxModel("One Factor", type = "RAM", 
	manifestVars = manifests, latentVars = latents, 
	mxPath(from = latents, to = manifests),
	mxPath(from = manifests, arrows = 2, values=.2),
	mxPath(from = latents, arrows = 2, free = FALSE, values = 1.0),
	mxPath(from = "one", to = manifests),
	mxData(cov(demoOneFactor), type='cov',
	       numObs=nrow(demoOneFactor), means=colMeans(demoOneFactor))
)

fm <- omxSetParameters(m1, "One Factor.M[1,1]", values = .1, free=FALSE)

expect_error(mxPowerSearch(m1, fm),
             "contains 'cov' data")
