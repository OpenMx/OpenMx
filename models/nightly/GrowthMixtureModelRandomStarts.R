#
#   Copyright 2007-2012 The OpenMx Project
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

# -----------------------------------------------------------------------
# Program: GrowthMixtureModelRandomStarts.R  
# Author: Unknown
# Date: 9999.99.99 
#
# ModelType: Growth Mixture
# DataType: Continuous
# Field: None
#
# Purpose: 
#      Growth mixture model with random starts
# 
#
# RevisionHistory:
#      Ross Gore -- 2011.06.16 updated & reformatted
# -----------------------------------------------------------------------


demo(GrowthMixtureModel_MatrixRaw)




trials <- 10
# how many trials?
# -------------------------------------


parNames <- names(omxGetParameters(gmm))
# place all of the parameter 
# names in a vector
# -------------------------------------


input <- matrix(NA, trials, length(parNames))
dimnames(input) <- list(c(1: trials), c(parNames))

output <- matrix(NA, trials, length(parNames))
dimnames(output) <- list(c(1: trials), c(parNames))

fit <- matrix(NA, trials, 5)
dimnames(fit) <- list(c(1: trials), c("Minus2LL", "Status", "Iterations", "pclass1", "time"))
# make matrices to hold everything
# -------------------------------------

input[,"p1"] <- runif(trials, 0.1, 0.9)
input[,"p1"] <- input[,"p1"]/(1-input[,"p1"])
# populate the class probabilities
# -------------------------------------


v <- c("vari1", "vars1", "vari2", "vars2", "residual")
input[,v] <- runif(trials*5, 0, 10)
# populate the variances
# -------------------------------------


m <- c("meani1", "means1", "meani2", "means2")
input[,m] <- runif(trials*4, -5, 5)
# populate the means
# -------------------------------------


r <- runif(trials*2, -0.9, 0.9)
scale <- c(
    sqrt(input[,"vari1"]*input[,"vars1"]),
    sqrt(input[,"vari2"]*input[,"vars2"]))
input[,c("cov1", "cov2")] <- r * scale
# populate the covariances
# -------------------------------------

for (i in 1: trials){
	temp1 <- omxSetParameters(gmm,
		labels=parNames,
		values=input[i,]
		)
		
	temp1@name <- paste("Starting Values Set", i)
		
	temp2 <- mxRun(temp1, unsafe=TRUE, suppressWarnings=TRUE)
	
	output[i,] <- omxGetParameters(temp2)
	fit[i,] <- c(
		temp2@output$Minus2LogLikelihood,
		temp2@output$status[[1]],
		temp2@output$iterations,
		round(temp2$classProbs@result[1,1], 4),
		temp2@output$wallTime
		)
	}
	
fit
table(round(fit[,1], 3), fit[,2])

# Serial Respecification
# -----------------------------------------------------------------------------




#require(snowfall)
#sfInit(parallel=TRUE, cpus=4)
#sfLibrary(OpenMx)

topModel <- mxModel("Top")	

makeModel <- function(modelNumber){
	temp <- mxModel(gmm, 
		independent=TRUE,
		name=paste("Iteration", modelNumber, sep=""))
	temp <- omxSetParameters(temp,
		labels=parNames,
		values=input[modelNumber,])
	return(temp)
}
	
mySubs <- lapply(1:trials, makeModel)
	
topModel@submodels <- mySubs

results <- mxRun(topModel, suppressWarnings=TRUE)

fitStats <- function(model){
	retval <- c(
		model@output$Minus2LogLikelihood,
		model@output$status[[1]],
		model@output$iterations,
		round(model$classProbs@result[1,1], 4)
		)	
	return(retval)
}

resultsFit <- t(omxSapply(results@submodels, fitStats))
#sfStop()

# Parallel Respecification
# -----------------------------------------------------------------------------


sum(fit[,5])
# Serial Time, in seconds (ignoring the overhead caused by the for loop)
# -------------------------------------

results@output$wallTime
# Parallel Time, in seconds
# -------------------------------------

# Compare Model Estimation Times
# -----------------------------------------------------------------------------


library(OpenMx)

# Reload OpenMx to offload the snowfall OpenMx library
# -----------------------------------------------------------------------------