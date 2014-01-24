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
# Program: BootstrapParallel.R  
# Author: Unknown
#  Date: 9999.99.99 
#
# ModelType: Parallel
# DataType: Continuous
# Field: None
#
# Purpose:
#      Bootstrap parallel models
#
# RevisionHistory:
#      Ross Gore -- 2011.06.16 updated & reformatted
# -----------------------------------------------------------------------

require(OpenMx)


lambda <- matrix(c(.8, .5, .7, 0), 4, 1)
nObs <- 500
nReps <- 10
nVar <- nrow(lambda)
specifics <- diag(nVar)
chl <- chol(lambda %*% t(lambda) + specifics)
# parameters for the simulation: lambda = factor loadings,
# specifics = specific variances
# -----------------------------------------------------------------------------


pStrt <- 3
pEnd <- pStrt + 2*nVar - 1
hStrt <- pEnd + 1
hEnd <- hStrt + 2*nVar - 1
# indices for parameters and hessian estimate in results
# -----------------------------------------------------------------------------


dn <- list()
dn[[1]] <- paste("Var", 1:4, sep="")
dn[[2]] <- dn[[1]]
# dimension names for OpenMx
# -----------------------------------------------------------------------------


randomCov <- function(nObs, nVar, chl, dn) {
  x <- matrix(rnorm(nObs*nVar), nObs, nVar)
  x <- x %*% chl
  thisCov <- cov(x)
  dimnames(thisCov) <- dn
  return(thisCov)  
}
# function to get a covariance matrix
# -----------------------------------------------------------------------------

createNewModel <- function(index, prefix, model) {
	modelname <- paste(prefix, index, sep='')
	data <- mxData(randomCov(nObs, nVar, chl, dn), type="cov", numObs=nObs)
	model <- mxModel(model, data)
	model <- mxRename(model, modelname)
	return(model)
}

getStats <- function(model) {
	retval <- c(model@output$status[[1]],
		max(abs(model@output$gradient)),
		model@output$estimate,
		sqrt(diag(solve(model@output$hessian))))
	return(retval)
}



obsCov <- randomCov(nObs, nVar, chl, dn)
# initialize obsCov for MxModel
# -----------------------------------------------------------------------------

results <- matrix(0, nReps, hEnd)
dnr <- c("inform", "maxAbsG", paste("lambda", 1:nVar, sep=""),
         paste("specifics", 1:nVar, sep=""),
         paste("hessLambda", 1:nVar, sep=""),
         paste("hessSpecifics", 1:nVar, sep=""))
dimnames(results)[[2]] <- dnr
# results matrix: get results for each simulation
# -----------------------------------------------------------------------------


template <- mxModel("stErrSim",
                       mxMatrix(name="lambda", type="Full", nrow=4, ncol=1,
                                free=TRUE, values=c(.8, .5, .7, 0)),
                       mxMatrix(name="specifics", type="Diag", nrow=4,
                                free=TRUE, values=rep(1, 4)),
                       mxAlgebra(lambda %*% t(lambda) + specifics,
                                 name="preCov", dimnames=dn),
                       mxData(observed=obsCov, type="cov", numObs=nObs),
                       mxMLObjective(covariance='preCov'),
                       independent = TRUE)
# instantiate MxModel
# -----------------------------------------------------------------------------

topModel <- mxModel("container")

submodels <- lapply(1:nReps, createNewModel, 'stErrSim', template)

names(submodels) <- imxExtractNames(submodels)
topModel@submodels <- submodels

modelResults <- mxRun(topModel, silent=TRUE, suppressWarnings=TRUE)

results <- t(omxSapply(modelResults@submodels, getStats))


results2 <- data.frame(results[which(results[,1] <= 1),])
# get rid of bad covergence results
# -----------------------------------------------------------------------------

means <- colMeans(results2)
stdevs <- sapply(results2, sd)
sumResults <- data.frame(matrix(dnr[pStrt:pEnd], 2*nVar, 1,
                                dimnames=list(NULL, "Parameter")))
sumResults$mean <- means[pStrt:pEnd]
sumResults$obsStDev <- stdevs[pStrt:pEnd]
sumResults$meanHessEst <- means[hStrt:hEnd]
sumResults$sqrt2meanHessEst <- sqrt(2) * sumResults$meanHessEst
# summarize the results
# -----------------------------------------------------------------------------

print(sumResults)
# print results
# -----------------------------------------------------------------------------

