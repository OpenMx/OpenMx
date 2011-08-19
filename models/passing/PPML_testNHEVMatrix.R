## Non-homogeneous Error Variance Test Model 2, Matrix Specified
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

factorModel <- mxModel("Two Factor NHEV, Matrix Specification",
	  mxData(dataCov, type="cov", numObs=100),
	  mxMatrix(
		type='Full',
		nrow=5,
		ncol=5,
		free=FALSE,
		values=c(0,0,0,2,1,
				 0,0,0,1,1,
				 0,0,0,1,1,
				 0,0,0,0,0,
				 0,0,0,0,0),
		byrow=TRUE,
		name='A'),
	
	  mxMatrix(
		type='Full',
		nrow=5,
		ncol=5,
		free=c(TRUE, TRUE, FALSE, FALSE, FALSE,
		       TRUE, TRUE, FALSE, FALSE, FALSE,
			   FALSE, FALSE, TRUE, FALSE, FALSE,
			   FALSE, FALSE, FALSE, TRUE, FALSE,
			   FALSE, FALSE, FALSE, FALSE, TRUE),
		values=c(1.0, 0.1, 0, 0, 0,
				 0.1, 2.0, 0, 0, 0,
				 0.0, 0, 3.0, 0, 0,
				 0, 0, 0, 1, 0,
				 0, 0, 0, 0, 1),
		labels=c("E1", "E4", NA, NA, NA,
				 "E4", "E2", NA, NA, NA,
				 NA, NA, "E3", NA, NA,
				 NA, NA, NA, NA, NA,
				 NA, NA, NA, NA, NA),
		byrow=TRUE,
		name='S'),
	mxMatrix(
		type='Full',
		nrow=3,
		ncol=5,
		free=FALSE,
		values=c(1, 0, 0, 0, 0,
				 0, 1, 0, 0, 0,
				 0, 0, 1, 0, 0),
		byrow=TRUE,
		name="F"),
	mxConstraint('E2' == 2*'E1', name="C1"),
	mxConstraint('E3' == 3*'E1', name="C2"),
	mxConstraint('E4' == 0.1*'E1', name="C3"),
	mxRAMObjective("A","S","F",dimnames=c("X", "Y", "Z", "G", "H"))
)

# Get results from un-PPMLed model
res1 <- mxRun(factorModel, suppressWarnings = TRUE)
# Get results from PPMLed model
res2 <- mxRun(imxTransformModelPPML(factorModel), suppressWarnings = TRUE)
# Make sure transform was applied
omxCheckTrue( (length(grep("(PPML Transformed)", res2@name, fixed = TRUE)) > 0) )
omxCheckTrue( is(res2@objective, "MxAlgebraObjective") )
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
