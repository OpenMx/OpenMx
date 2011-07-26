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
require(MASS)
dataCov <- rbind(c(4,2,2),c(2,4,2),c(2,2,4))
manifests <- c('X','Y','Z')
colnames(dataCov) <- manifests
rownames(dataCov) <- manifests
dataTest <- mvrnorm(n=100, c(0,0,0), dataCov)
latents <- c('G', 'H')
factorModel <- mxModel("Two Factor NHEV",
      type="RAM",
      manifestVars = manifests,
      latentVars = latents,
      
	  mxPath(from='G', to=manifests,value=c(2,1,1),free=FALSE),
	  mxPath(from='H', to=manifests,value=c(1,1,1),free=FALSE),
	  
	  mxConstraint('E2' == 0.1*'E1', name="Constraint"),
      mxPath(from=manifests, arrows=2,value=c(10,10,10), labels=c('E1','E1','E1')),	  
	  mxPath(from='X', to='Y', arrows=2, value=1, labels='E2'),
	  
      mxPath(from=latents, arrows=2,
            values=1.0),
	  mxPath(from="one", to=manifests, arrows=1, values=0, free=FALSE),
      mxData(dataTest, type="raw", numObs=100)
			)

# Get results from un-PPMLed model
res1 <- mxRun(factorModel)
# Get results from PPMLed model
res2 <- mxRun(imxTransformModelPPML(factorModel))
# Make sure transform was applied
omxCheckTrue( (length(grep("(PPML Transformed)", res2@name, fixed = TRUE)) > 0) )
omxCheckTrue( is(res2@objective, "MxAlgebraObjective") )
# Check error param, latent variances
omxCheckCloseEnough(res2@output$estimate['_PPML_NHEV_ErrParam'], res1@output$estimate['E1'], .001)
omxCheckCloseEnough(res2@output$estimate[2], res1@output$estimate[3], .001)
omxCheckCloseEnough(res2@output$estimate[3], res1@output$estimate[4], .001)
# Get restored results
res3 <- imxRestoreResultPPML(factorModel, res1)
# Check log likelihood, error variances, latent variances
omxCheckCloseEnough(res3@output$minimum, res1@output$minimum, .001)
omxCheckCloseEnough(res3@output$estimate['E1'], res1@output$estimate['E1'], .001)
omxCheckCloseEnough(res3@output$estimate['E2'], res1@output$estimate['E2'], .001)
omxCheckCloseEnough(res3@output$estimate[3], res1@output$estimate[3], .001)
omxCheckCloseEnough(res3@output$estimate[4], res1@output$estimate[4], .001)
