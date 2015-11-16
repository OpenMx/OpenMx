## Non-homogeneous Error Variance Test Model
# Two factor model
# Two latents predict three manifests.  Regression weights between latents and
# manifests are fixed, variances for latents are free, variances between manifests
# are specified by a fixed error covariance matrix multiplied by a scalar error
# parameter.  This is implemented by labeling all free terms in the symmetric matrix
# and linearly constraining them relative to each other.
# In this model, the error matrix is very simple: homogeneous along the diagonal, but
# with a covariance between two of the manifest variables.
# Written by Daniel Hackett, July 2011
require(OpenMx)


manifests <- c('X','Y','Z')
latents <- c('G', 'H')

dataCov <- rbind(c(4,3,1),c(3,4,1),c(1,1,2))
dataMean <- c(2,1,3)
#set.seed(42)
#raw <- mvtnorm::rmvnorm(n=100, mu=dataMean, Sigma=dataCov)
names(dataMean) <- manifests
#colnames(raw) <- manifests
colnames(dataCov) <- manifests
rownames(dataCov) <- manifests

factorModel <- mxModel("Two Factor NHEV",
      type="RAM",
      manifestVars = manifests,
      latentVars = latents,
      
	  mxPath(from='G', to=manifests,value=c(3,2,1),free=FALSE),
	  mxPath(from='H', to=manifests,value=c(1,2,3),free=FALSE),

	  mxConstraint('E2' == 0.5*'E1', name="Constraint"),
      mxPath(from=manifests, arrows=2,value=c(10,10,10), labels=c('E1','E1','E1')),	  
	  mxPath(from='X', to='Y', arrows=2, value=5, labels='E2'),
		mxPath(from='one', to=latents, arrows=1, value=c(1,1), labels=c('m_G', 'm_H')),

		mxPath(from='G', to='H', arrows=2, values=0.1, labels='C_GH', free=TRUE),	  
      mxPath(from=latents, arrows=2,
            values=1.0),
	  mxData(dataCov, type="cov", means=dataMean, numObs=100)
#		mxData(raw, type="raw", numObs=100)
			)
imxPPML.Test.Battery(factorModel, testMissingness = FALSE, tolerances = c(NA, .002, .001) )
