## Third simplest possible PPML candidate model
# Two factor model
# Two latents predict three manifests.  Regression weights between latents and
# manifests are fixed, variances for all five variables are free.
# Written by Daniel Hackett, June 2011

require(OpenMx)
dataTest <- rbind(c(4,2,2),c(2,4,2),c(2,2,4))
manifests <- c('X','Y','Z')
colnames(dataTest) <- manifests
rownames(dataTest) <- manifests
latents <- c('G', 'H')
factorModel <- mxModel("Two Factor",
      type="RAM",
      manifestVars = manifests,
      latentVars = latents,
      mxPath(from='G', to=manifests,value=c(2,1,1),free=FALSE),
	  mxPath(from='H', to=manifests,value=c(1,1,1),free=FALSE),
      mxPath(from=manifests, arrows=2,value=c(1,1,1), labels=c('E1','E1','E1')),
      mxPath(from=latents, arrows=2,
            values=1.0),
      mxData(dataTest, type="cov",
            numObs=100))

imxTestPPML(factorModel)
