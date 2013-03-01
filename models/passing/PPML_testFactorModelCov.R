## Third simplest possible PPML candidate model
# Two factor model
# Two latents predict three manifests.  Regression weights between latents and
# manifests are fixed, variances for all five variables are free.
# Written by Daniel Hackett, June 2011

require(OpenMx)
dataTest <- rbind(c(4,2,2),c(2,4,2),c(2,2,4))
dataMean <- c(2, 1, 0)

manifests <- c('X','Y','Z')
colnames(dataTest) <- manifests
rownames(dataTest) <- manifests
names(dataMean) <- manifests

latents <- c('G', 'H')
factorModel <- mxModel("Two Factor",
      type="RAM",
      manifestVars = manifests,
      latentVars = latents,
	  
	  # Regression loadings
      mxPath(from='G', to=manifests,value=c(3,2,1),free=FALSE),
	  mxPath(from='H', to=manifests,value=c(1,1,1),free=FALSE),
	  
	  # Means
	  mxPath(from="one", to=latents, value=c(0,0), free=TRUE),
	  
	  # Covariances
	  mxPath(from=manifests, arrows=2,value=c(1,1,1), labels=c('E1','E1','E1')),
      mxPath(from=latents, arrows=2,value=c(0,2)),
	  mxPath(from='G', to='H', arrows=2, value=0, free=TRUE), # Saturate
	  
      mxData(dataTest, type="cov", means=dataMean,
            numObs=100))

imxTestPPML(factorModel)
