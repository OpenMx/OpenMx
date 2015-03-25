require(OpenMx)
require(MASS)

# Seed
set.seed(4)

# Data covariance matrix
dataCov <- rbind(c(10,3,1),
				 c(3,8,2),
				 c(1,2,6))
dataMean <- c(1,1,10)
#dataMean <- c(2,-2,1)

# Variable names
manifests <- c('X','Y','Z')
latents <- c('G', 'H')

# Dimnames
colnames(dataCov) <- manifests
rownames(dataCov) <- manifests
names(dataMean) <- manifests

factorModel <- mxModel("Factor Model 2L3M",
      type="RAM",
      manifestVars = manifests,
      latentVars = latents,
	  
	  # Regression Loadings
      mxPath(from='G', to=manifests,value=c(1,2,3),free=FALSE),
	  mxPath(from='H', to=manifests,value=c(1,1,1),free=FALSE),
	  
	  # Variances
      mxPath(from=manifests, arrows=2,value=c(1,1,1), labels=c('E1','E1','E1')), # error variance
      mxPath(from=latents, arrows=2,values=1.0, labels=c("VG", "VH")), # variances of the latents
	 
	  # Latent means vector
      mxPath(from="one", to=latents, arrows=1, values=0, labels=c("MG", "MH"), free=TRUE),

	  mxData(dataCov, mean=dataMean, type="cov",numObs=5)
	  )

# Latent Covariance	Path
satPath <- mxPath(from='G', to='H', labels="Cov", arrows=2, values=0.1)	# Saturated (covariance between the latents)

for (fixALatent in 0:1)
{
	testModel <- factorModel
	if (as.logical(fixALatent))
		testModel <- mxModel(testModel, mxPath(from="one", to=latents, arrows=1, values=0, labels=c("MG", "MH"), free=c(TRUE, FALSE)) )

	# Test unsaturated latent covariance matrix
	imxPPML.Test.Battery(testModel, testPermutations=FALSE, verbose=TRUE, tolerances=c(NA, .0001, .0001) ) # NA -> Don't check covariance data w/ means
	
	# Test saturated latent covariance matrix
	testModel <- mxModel(testModel, satPath)
	imxPPML.Test.Battery(testModel, testPermutations=FALSE, verbose=TRUE, tolerances=c(NA, .0001, .0001) ) # NA -> Don't check covariance data w/ means
}
