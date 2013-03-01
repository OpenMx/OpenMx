## Third simplest possible PPML candidate model
# Two factor model
# Two latents predict three manifests.  Regression weights between latents and
# manifests are fixed, variances for all five variables are free.
# Written by Daniel Hackett, June 2011

require(OpenMx)
require(MASS)
set.seed(4)
#dataCov <- rbind(c(3,3,4),c(3,6,7),c(4,7,11))
dataCov <- rbind(c(4,3,2,1), c(3,4,3,2), c(2,3,4,3), c(1,2,3,4))
manifests <- c('X','Y','Z','M')
colnames(dataCov) <- manifests
rownames(dataCov) <- manifests
dataTest <- mvrnorm(n=5000, c(2,-2,1,-1), dataCov)
latents <- c('G', 'H')
factorModel <- mxModel("Two Factor",
      type="RAM",
      manifestVars = manifests,
      latentVars = latents,
	  
      mxPath(from='G', to=manifests,value=c(1,2,3),free=FALSE),
	  mxPath(from='H', to=manifests,value=c(1,1,1),free=FALSE),
	  
      mxPath(from=manifests, arrows=2,value=c(1,1,1), labels=c('E1','E1','E1')), # error variance
      mxPath(from=latents, arrows=2,values=1.0), # variances of the latents
	  
	  mxPath(from='G', to='H', arrows=2,values=0.1),	# Saturated (covariance between the latents)
	 
      mxPath(from="one", to=latents, arrows=1, values=0, free=TRUE),
      mxData(dataTest, type="raw",
            numObs=5000))

imxTestPPML(factorModel)
