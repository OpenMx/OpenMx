## Non-homogeneous Error Variance Test Model 2
# Two factor model
# Two latents predict three manifests.  Regression weights between latents and
# manifests are fixed, variances for latents are free, variances between manifests
# are specified by a fixed error covariance matrix multiplied by a scalar error
# parameter.  This is implemented by labeling all free terms in the symmetric matrix
# and linearly constraining them relative to each other.
# In this case, unlike in the simpler NHEV test model, this error matrix is not
# homogeneous on the diagonal.
# Written by Daniel Hackett, July 2011

require(OpenMx)
require(MASS)
dataCov <- rbind(c(4,2,2),c(2,4,2),c(2,2,4))
manifests <- c('X','Y','Z')
colnames(dataCov) <- manifests
rownames(dataCov) <- manifests
latents <- c('G', 'H')
factorModel <- mxModel("Two Factor NHEV",
      type="RAM",
      manifestVars = manifests,
      latentVars = latents,
      mxPath(from='G', to=manifests,value=c(2,1,1),free=FALSE),
	  mxPath(from='H', to=manifests,value=c(1,1,1),free=FALSE),
	  
	  mxPath(from='X', arrows=2, value=1, labels='E1'),
	  mxPath(from='Y', arrows=2, value=2, labels='E2'),
	  mxPath(from='Z', arrows=2, value=3, labels='E3'),
	  mxPath(from='X', to='Y', arrows=2, value =0.1, labels='E4'),
	  mxConstraint('E2' == 2*'E1', name="C1"),
	  mxConstraint('E3' == 3*'E1', name="C2"),
	  mxConstraint('E4' == 0.1*'E1', name="C3"),
      mxPath(from=latents, arrows=2,
            values=1.0),
	  mxData(dataCov, type="cov", numObs=100)
			)

# Get results from un-PPMLed model
res1 <- mxRun(factorModel)
# Get results from PPMLed model
res2 <- mxRun(imxTransformModelPPML(factorModel))
# Check error param, latent variances
omxCheckCloseEnough(res2@output$estimate['_PPML_NHEV_ErrParam'], res1@output$estimate['E1'], .001)
omxCheckCloseEnough(res2@output$estimate[2], res1@output$estimate[5], .001)
omxCheckCloseEnough(res2@output$estimate[3], res1@output$estimate[6], .001)
# Get restored results
res3 <- imxRestoreResultPPML(factorModel, res1)
# Check log likelihood, error variances, latent variances
omxCheckCloseEnough(res3@output$minimum, res1@output$minimum, .001)
omxCheckCloseEnough(res3@output$estimate['E1'], res1@output$estimate['E1'], .001)
omxCheckCloseEnough(res3@output$estimate['E2'], res1@output$estimate['E2'], .001)
omxCheckCloseEnough(res3@output$estimate['E3'], res1@output$estimate['E3'], .001)
omxCheckCloseEnough(res3@output$estimate['E4'], res1@output$estimate['E4'], .001)
omxCheckCloseEnough(res3@output$estimate[5], res1@output$estimate[5], .001)
omxCheckCloseEnough(res3@output$estimate[6], res1@output$estimate[6], .001)
