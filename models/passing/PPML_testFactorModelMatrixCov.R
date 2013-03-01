## Two Factor Model -- Matrix Specification
# One factor model
# One latent predicts three manifests.  Regression weights between latents and
# manifests are fixed, variances for latents are free, error variances of the
# manifest variables are homogeneous.
# Written by Daniel Hackett, August 2011

require(OpenMx)
dataCov <- rbind(c(4,2,2),c(2,4,2),c(2,2,4))
manifests <- c('X','Y','Z')
colnames(dataCov) <- manifests
rownames(dataCov) <- manifests

latents <- c('G', 'H')

factorModel <- mxModel("Two Factor Matrix",
	  mxMatrix(
		type = "Full",
		nrow = 5,
		ncol = 5,
		free = FALSE,
		values = c(0, 0, 0, 3, 4,
				   0, 0, 0, 2, 3,
				   0, 0, 0, 1, 2,
				   0, 0, 0, 0, 0,
				   0, 0, 0, 0, 0),
		byrow = TRUE,
		name = "A"),
	  mxMatrix(
		type = "Symm",
		nrow = 5,
		ncol = 5,
		free = c(T, F, F, F, F,
				 F, T, F, F, F,
				 F, F, T, F, F,
				 F, F, F, T, T,
				 F, F, F, T, T),
		values = c(1, 0, 0, 0, 0,
				   0, 1, 0, 0, 0,
				   0, 0, 1, 0, 0,
				   0, 0, 0, 1, 0.1,
				   0, 0, 0, 0.1, 1),
		labels = c("residual", NA, NA, NA, NA,
					NA, "residual", NA, NA, NA,
					NA, NA, "residual", NA, NA,
					NA, NA, NA, "VarG", "CovGH",
					NA, NA, NA, "CovGH", "VarH"),
		byrow = TRUE,
		name = "S"),
	  mxMatrix(
		type = "Full",
		nrow = 3,
		ncol = 5,
		free = FALSE,
		values = c(1, 0, 0, 0, 0,
				   0, 1, 0, 0, 0,
				   0, 0, 1, 0, 0),
		byrow = TRUE,
		name="F"),
		
      mxData(dataCov, type="cov", numObs=140),
	  mxRAMObjective("A","S","F",dimnames=c("X", "Y", "Z", "G","H"))
			)
imxTestPPML(factorModel)
