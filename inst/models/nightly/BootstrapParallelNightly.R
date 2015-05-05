#
#   Copyright 2007-2010 The OpenMx Project
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
require(snowfall)
sfInit(parallel = TRUE, cpus = 8)
sfLibrary(OpenMx)

set.seed(10)

# parameters for the simulation: lambda = factor loadings,
# specifics = specific variances
lambda <- matrix(c(.8, .5, .7, 0), 4, 1)
nObs <- 500
nReps <- 1000
nVar <- nrow(lambda)
specifics <- diag(nVar)
chl <- chol(lambda %*% t(lambda) + specifics)

# indices for parameters and hessian estimate in results
pStrt <- 3
pEnd <- pStrt + 2*nVar - 1
hStrt <- pEnd + 1
hEnd <- hStrt + 2*nVar - 1

# dimension names for OpenMx
dn <- list()
dn[[1]] <- paste("Var", 1:4, sep="")
dn[[2]] <- dn[[1]]

# function to get a covariance matrix
randomCov <- function(nObs, nVar, chl, dn) {
  x <- matrix(rnorm(nObs*nVar), nObs, nVar)
  x <- x %*% chl
  thisCov <- cov(x)
  dimnames(thisCov) <- dn
  return(thisCov)  
}

createNewModel <- function(model, modelname) {
	model$data$observed <- randomCov(nObs, nVar, chl, dn)
	model <- mxModel(name=modelname, model)
	return(model)
}

getStats <- function(model) {
	retval <- c(code=model$output$status[[1]],
		grad=sqrt(sum(model$output$gradient^2)),
		model$output$estimate,
		sqrt(2*diag(solve(model$output$hessian))))
	return(retval)
}


# initialize obsCov for MxModel
obsCov <- randomCov(nObs, nVar, chl, dn)

# results matrix: get results for each simulation
results <- matrix(0, nReps, hEnd)
dnr <- c("inform", "normG", paste("lambda", 1:nVar, sep=""),
         paste("specifics", 1:nVar, sep=""),
         paste("hessLambda", 1:nVar, sep=""),
         paste("hessSpecifics", 1:nVar, sep=""))
dimnames(results)[[2]] <- dnr

# instantiate MxModel
template <- mxModel(name="stErrSim",
                       mxMatrix(name="lambda", type="Full", nrow=4, ncol=1,
                                free=TRUE, values=c(.8, .5, .7, 0)),
                       mxMatrix(name="specifics", type="Diag", nrow=4,
                                free=TRUE, values=rep(1, 4)),
                       mxAlgebra(lambda %*% t(lambda) + specifics,
                                   name="preCov", dimnames=dn),
                         mxData(observed=obsCov, type="cov", numObs=nObs),
                      mxExpectationNormal(covariance='preCov'),
                    mxFitFunctionML(),
                    independent = TRUE)

submodels <- lapply(1:nReps, function(x) {
  createNewModel(template, paste('stErrSim', x, sep=''))
})

names(submodels) <- imxExtractNames(submodels)

topModel <- mxModel('container', submodels)

modelResults <- mxRun(topModel, silent=TRUE, suppressWarnings=TRUE)

results <- t(omxSapply(modelResults$submodels, getStats))

sfStop()

# get rid of bad covergence results
results2 <- results[which(results[,"code"] <= 1),]

# summarize the results
means <- colMeans(results2)
stdevs <- apply(results2, 2, sd)
sumResults <- data.frame(matrix(dnr[pStrt:pEnd], 2*nVar, 1,
                                dimnames=list(NULL, "Parameter")))
sumResults$mean <- means[pStrt:pEnd]
sumResults$obsStDev <- stdevs[pStrt:pEnd]
sumResults$meanHessEst <- means[hStrt:hEnd]

print(sumResults)

omxCheckCloseEnough(means["grad"], 0, .05)
omxCheckCloseEnough(sumResults$mean, c(lambda, diag(specifics)), .05)
omxCheckCloseEnough(sumResults$obsStDev, sumResults$meanHessEst, .05)

detach("package:snowfall")
unloadNamespace("snowfall")
detach("package:snow")
