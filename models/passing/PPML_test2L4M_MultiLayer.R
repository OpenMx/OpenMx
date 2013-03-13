require(OpenMx)
require(MASS)

# Seed
set.seed(4)

# Data covariance matrix
dataCov <- rbind(c(4,3,2,1), 
				 c(3,4,3,2),
				 c(2,3,4,3),
				 c(1,2,3,4))
dataMean <- c(2,-2,1,-1)

# Variable names
manifests <- c('X','Y','Z','W')
latents <- c('G', 'H', 'root')

# Dimnames
colnames(dataCov) <- manifests
rownames(dataCov) <- manifests
names(dataMean) <- manifests

factorModel <- mxModel("Factor Model 2L4M",
      type="RAM",
      manifestVars = manifests,
      latentVars = latents,
	  
	  # Regression Loadings
			mxPath(from='root', to=c('G', 'H'), value=c(-1,1), free=FALSE),
      mxPath(from='G', to=manifests,value=c(1,2,3,4),free=FALSE),
		  mxPath(from='H', to=manifests,value=c(1,1,1,1),free=FALSE),
	  
	  # Variances
      mxPath(from=manifests, arrows=2,value=c(1,1,1,1), labels=c('E1','E1','E1','E1')), # error variance
      mxPath(from=latents, arrows=2,values=1.0, labels=c("VG", "VH", "Vroot")), # variances of the latents
	 
	  # Latent means vector
      mxPath(from="one", to=latents, arrows=1, values=0, labels=c("MG", "MH","Mroot"), free=TRUE),

	  mxData(dataCov, mean=dataMean, type="cov",numObs=100)
	  )

# Latent Covariance	Path
satPath <- mxPath(from='G', to='H',labels="Cov", arrows=2,values=0.1)	# Saturated (covariance between the latents)

for (fixALatent in 0:1)
{
	testModel <- factorModel
	if (as.logical(fixALatent))
		testModel <- mxModel(testModel, mxPath(from="one", to=latents, arrows=1, values=0, labels=c("MG", "MH"), free=c(TRUE, FALSE)) )
	
	# Test unsaturated latent covariance matrix
	imxPPML.Test.Battery(testModel, testPermutations=FALSE, testMissingness=FALSE, testEstimates=FALSE)
	gc()
	# Test saturated latent covariance matrix
	testModel <- mxModel(testModel, satPath)
	imxPPML.Test.Battery(testModel, testPermutations=FALSE, testMissingness=FALSE, testEstimates=FALSE)
	gc()
}
