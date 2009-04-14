#Latent Growth Curve - Matrix Specification
#Needs mean structure components
#Returning Error on mxRun: details below (R. Estabrook, 09Feb26)

require(OpenMx)

# Data
wisc.cov<-matrix(c(33.73019, 25.46499, 30.88961, 40.51781,
	25.46499, 37.28895, 33.82058, 47.40553,
	30.88961, 33.82058, 53.57812, 62.25334,
	40.51781, 47.40553, 62.25334, 113.74121),
	nrow=4, ncol=4, byrow=TRUE)

# Define the model
model <- mxModel()
model <- mxModel(model, mxMatrix("Full", name = "A", nrow = 6, ncol = 6))
model <- mxModel(model, mxMatrix("Symm", name = "S", nrow = 6, ncol = 6))
model <- mxModel(model, mxMatrix("Full", name = "F", nrow = 4, ncol = 6))

# Specify "A" Matrix
model[["A"]]@spec[6,2] <- "Basis Loading V2" #Latent Basis Slope Loading (Free)
model[["A"]]@spec[6,3] <- "Basis Loading V4" #Latent Basis Slope Loading (Free)

# Values for "A" Matrix
model[["A"]]@values[5,1] <- 1 #Intercept Loading
model[["A"]]@values[5,2] <- 1 #Intercept Loading
model[["A"]]@values[5,3] <- 1 #Intercept Loading
model[["A"]]@values[5,4] <- 1 #Intercept Loading
model[["A"]]@values[6,4] <- 1 #Latent Basis Slope Loading (Fixed)

model[["A"]]@values[6,2] <- .2 #Latent Basis Slope Loading (Free)
model[["A"]]@values[6,3] <- .6 #Latent Basis Slope Loading (Free)

# Specify "S" Matrix
model[["S"]]@spec[1,1] <- "Manifest Residual"
model[["S"]]@spec[2,2] <- "Manifest Residual"
model[["S"]]@spec[3,3] <- "Manifest Residual"
model[["S"]]@spec[4,4] <- "Manifest Residual"
model[["S"]]@spec[5,5] <- "Latent Intercept Variance"
model[["S"]]@spec[6,6] <- "Latent Slope Variance"
model[["S"]]@spec[5,6] <- "Latent Covariance"

# Values for "S" Matrix
model[["S"]]@values[1,1] <- 11
model[["S"]]@values[2,2] <- 11
model[["S"]]@values[3,3] <- 11
model[["S"]]@values[4,4] <- 11
model[["S"]]@values[5,5] <- 20
model[["S"]]@values[6,6] <- 24
model[["S"]]@values[5,6] <- 14

# Specify "F" Matrix

# Values for "F" Matrix
model[["F"]]@values[1,1] <- 1
model[["F"]]@values[2,2] <- 1
model[["F"]]@values[3,3] <- 1
model[["F"]]@values[4,4] <- 1

# Define the objective function
objective <- mxRAMObjective("A", "S", "F", "objective")

# Define the observed covariance matrix
covMatrix <- mxData(wisc.cov, type='cov', numObs = 100)

# Add the objective function and the data to the model
model <- mxModel(model, objective, covMatrix)

# Run the job
# Currently returns "Error in if (uni) "U" else "N" : missing value where TRUE/FALSE needed"
model <- mxRun(model)
print(model@output)

# Expected Results
	# Manifest Residual			11.005 (0.771)
	# Latent Intercept Variance	19.713 (0.394)
	# Latent Slope Variance		24.140 (0.575)
	# Latent Covariance			14.004 (3.050)
	# Basis Coefficient: V2		0.234 (0.012)
	# Basis Coefficient: V4		0.527 (0.011)
# Expected Fit Statistics (when mean structure is included)
	# Log Likelihood	-2491.070
	# AIC				4998.139
	# RMSEA				.116
