#
#   Copyright 2007-2009 The OpenMx Project
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

require(OpenMx)

#Latent Growth Curve - Matrix Specification
#Needs mean structure components
#Returning Error on mxRun: details below (R. Estabrook, 09Feb26)

# Data
wisc.cov<-matrix(c( 33.73019, 25.46499, 30.88961, 40.51781,
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
model[["A"]]@labels[2,6] <- "Basis Loading V2" #Latent Basis Slope Loading (Free)
model[["A"]]@labels[3,6] <- "Basis Loading V4" #Latent Basis Slope Loading (Free)
model[["A"]]@free[2,6] <- TRUE 
model[["A"]]@free[3,6] <- TRUE

# Values for "A" Matrix
model[["A"]]@values[1,5] <- 1 #Intercept Loading
model[["A"]]@values[2,5] <- 1 #Intercept Loading
model[["A"]]@values[3,5] <- 1 #Intercept Loading
model[["A"]]@values[4,5] <- 1 #Intercept Loading
model[["A"]]@values[4,6] <- 1 #Latent Basis Slope Loading (Fixed)

model[["A"]]@values[2,6] <- .2 #Latent Basis Slope Loading (Free)
model[["A"]]@values[3,6] <- .6 #Latent Basis Slope Loading (Free)

# Specify "S" Matrix
model[["S"]]@labels[1,1] <- "Manifest Residual"
model[["S"]]@labels[2,2] <- "Manifest Residual"
model[["S"]]@labels[3,3] <- "Manifest Residual"
model[["S"]]@labels[4,4] <- "Manifest Residual"
model[["S"]]@labels[5,5] <- "Latent Intercept Variance"
model[["S"]]@labels[6,6] <- "Latent Slope Variance"
model[["S"]]@labels[5,6] <- "Latent Covariance"
model[["S"]]@labels[6,5] <- "Latent Covariance"

model[["S"]]@free[1,1] <- TRUE
model[["S"]]@free[2,2] <- TRUE
model[["S"]]@free[3,3] <- TRUE
model[["S"]]@free[4,4] <- TRUE
model[["S"]]@free[5,5] <- TRUE
model[["S"]]@free[6,6] <- TRUE
model[["S"]]@free[5,6] <- TRUE
model[["S"]]@free[6,5] <- TRUE

# Values for "S" Matrix
model[["S"]]@values[1,1] <- 11
model[["S"]]@values[2,2] <- 11
model[["S"]]@values[3,3] <- 11
model[["S"]]@values[4,4] <- 11
model[["S"]]@values[5,5] <- 20
model[["S"]]@values[6,6] <- 24
model[["S"]]@values[5,6] <- 14
model[["S"]]@values[6,5] <- 14

# Values for "F" Matrix
model[["F"]]@values[1,1] <- 1
model[["F"]]@values[2,2] <- 1
model[["F"]]@values[3,3] <- 1
model[["F"]]@values[4,4] <- 1


# Define the objective function
objective <- mxRAMObjective("A", "S", "F")

# Define the observed covariance matrix
covMatrix <- mxData(wisc.cov, type='cov', numObs = 100)

# Add the objective function and the data to the model
model <- mxModel(model, objective, covMatrix)

# Run the job
model <- mxRun(model)

estimates <- model@output$estimate

# Results from old Mx:
omxCheckCloseEnough(estimates[['Basis Loading V2']], 0.1110, 0.01)
omxCheckCloseEnough(estimates[['Basis Loading V4']], 0.3343, 0.01)
omxCheckCloseEnough(estimates[['Manifest Residual']], 10.3184, 0.01)
omxCheckCloseEnough(estimates[['Latent Intercept Variance']], 24.2041, 0.01)
omxCheckCloseEnough(estimates[['Latent Covariance']], 16.8327, 0.01)
omxCheckCloseEnough(estimates[['Latent Slope Variance']], 45.8975, 0.01)

# MPlus Results?
# omxCheckCloseEnough(estimates[['Basis Loading V2']], 0.234, 0.01)
# omxCheckCloseEnough(estimates[['Basis Loading V4']], 0.527, 0.01)
# omxCheckCloseEnough(estimates[['Manifest Residual']], 11.005, 0.01)
# omxCheckCloseEnough(estimates[['Latent Intercept Variance']], 19.713, 0.01)
# omxCheckCloseEnough(estimates[['Latent Covariance']], 14.004, 0.01)
# omxCheckCloseEnough(estimates[['Latent Slope Variance']], 24.140, 0.01)

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
