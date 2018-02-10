#
#   Copyright 2007-2018 The OpenMx Project
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


#
#  OpenMx Mixed Ordinal & Continuous Data Example
#  Revision history:
#		Michael Neale 7 Jan 2011
#       Revised Timothy Brick 10 Jan 2011

require(OpenMx)


# First simulation: 2 by 2, one factor for each, uncorrelated.
# Because there's no correlation between the factors, the result of the joint optimization
# is the sum of the ordinal likelihood and the continuous likelihood.
# All estimates on both sides should match.

nOrdinalVariables <- 2
nContinuousVariables <- 2
nVariables <- nOrdinalVariables + nContinuousVariables
nFactors <- 2
nThresholds <- 1
nSubjects <- 200
useOptimizer <- FALSE

loadings <- matrix(c(.7,.7, 0, 0, 0, 0, .7, .7),nrow=nVariables,ncol=nFactors)
startLoads <- loadings * .5
startFree <- as.logical(startLoads * useOptimizer)
residuals <- 1 - diag(loadings %*% t(loadings))
sigma <- loadings %*% t(loadings) + vec2diag(residuals)

mu <- matrix(0,nrow=nVariables,ncol=1)

# Step 1: simulate multivariate normal data for first simulation
set.seed(1234)
continuousData <- data.frame(matrix(mvtnorm::rmvnorm(n=nSubjects,mu,sigma), nrow=nSubjects, ncol=nVariables))

# Step 2: chop continuous variables into ordinal data 
# with nThresholds+1 approximately equal categories, based on 1st variable
quants <- quantile(continuousData[,1],  probs = c((1:nThresholds)/(nThresholds+1)))
ordinalData <- matrix(0,nrow=nSubjects,ncol=nVariables)
for(i in 1:nVariables) {
   ordinalData[,i] <- cut(as.vector(continuousData[,i]),c(-Inf,quants,Inf))
}

# Step 3: make the ordinal variables into R factors and make a joint data frame with both variables in it, innit?
ordinalData <- mxFactor(as.data.frame(ordinalData),levels=c(1:(nThresholds+1)))
jointData <- data.frame(continuousData[,(nOrdinalVariables+1):nVariables],
	ordinalData[,1:nOrdinalVariables])

ordinalNames <- paste("IamOrdinal", 1:nOrdinalVariables, sep="")
continuousNames <- paste("IamContinuous", 1:nContinuousVariables, sep="")
jointNames <- c(continuousNames, ordinalNames)
minOrd <- 1 + nContinuousVariables
maxOrd <- nOrdinalVariables + nContinuousVariables
minCont <- 1
maxCont <- nContinuousVariables
ordCols <- 1:nOrdinalVariables + nContinuousVariables
contCols <- 1:nContinuousVariables
names(jointData) <- jointNames
names(continuousData) <- jointNames
names(ordinalData) <- jointNames 

checkPermutations <- function(model) {
	model$data$.needSort <- FALSE
	d1 <- model$data$observed
	result <- expand.grid(seed=1:10, jointConditionOn=c('ordinal', 'continuous'), fit=NA)
	for (rx in 1:nrow(result)) {
		#print(rx)
		set.seed(result[rx,'seed'])
		model$data$observed <- d1[sample.int(nrow(d1)),]
		rownames(model$data$observed) <- 1:nrow(d1)
		#print(model$data$observed[,model$expectation$dims])
		model$fitfunction$jointConditionOn <- as.character(result[rx,'jointConditionOn'])
		#model$fitfunction$verbose <- 3L
		got <- mxRun(mxModel(model, mxComputeOnce('fitfunction', 'fit')), silent=TRUE)
		result[rx,'fit'] <- got$output$fit
	}
	print(result)
	omxCheckEquals(length(table(round(result$fit,3))), 1)
}

# Step 4: Set up actual models, simulation 1

# Model spec for all continuous case

contModel <- mxModel("contModel", 
    mxMatrix("Full", nVariables, nFactors, values=startLoads, free=startFree, lbound=-.99, ubound=.99, name="L"),
    mxMatrix("Unit", nVariables, 1, name="vectorofOnes"),
    mxMatrix("Zero", 1, nVariables, name="M"),
    mxAlgebra(vectorofOnes - (diag2vec(L %*% t(L))) , name="E"),
    mxAlgebra(L %*% t(L) + vec2diag(E), name="impliedCovs"),
    mxMatrix("Full", 
           name="thresholdDeviations", nrow=nThresholds, ncol=nOrdinalVariables,
           values=.2,
           free = useOptimizer, 
           lbound = rep( c(-Inf,rep(.01,(nThresholds-1))) , nOrdinalVariables),
           dimnames = list(c(), ordinalNames)),
    mxMatrix("Lower",nThresholds,nThresholds,values=1,free=F,name="unitLower"),
    mxAlgebra(unitLower %*% thresholdDeviations, name="thresholdMatrix"),
    mxAlgebra(impliedCovs[minCont:maxCont,minCont:maxCont], name="useCov"),
    mxAlgebra(M[1,minCont:maxCont], name="useM"),
    mxFitFunctionML(),mxExpectationNormal(covariance="useCov", means="useM", dimnames = continuousNames),
    mxData(observed=jointData, type='raw')
)
checkPermutations(contModel)

# Model spec for all ordinal case

ordModel <- mxModel("ordModel", 
    mxMatrix("Full", nVariables, nFactors, values=startLoads, free=startFree, lbound=-.99, ubound=.99, name="L"),
    mxMatrix("Unit", nVariables, 1, name="vectorofOnes"),
    mxMatrix("Zero", 1, nVariables, name="M"),
    mxAlgebra(vectorofOnes - (diag2vec(L %*% t(L))) , name="E"),
    mxAlgebra(L %*% t(L) + vec2diag(E), name="impliedCovs"),
    mxMatrix("Full", 
           name="thresholdDeviations", nrow=nThresholds, ncol=nOrdinalVariables,
           values=.2,
           free = useOptimizer, 
           lbound = rep( c(-Inf,rep(.01,(nThresholds-1))) , nOrdinalVariables),
           dimnames = list(c(), ordinalNames)),
    mxMatrix("Lower",nThresholds,nThresholds,values=1,free=F,name="unitLower"),
    mxAlgebra(unitLower %*% thresholdDeviations, name="thresholdMatrix"),
    mxAlgebra(impliedCovs[minOrd:maxOrd,minOrd:maxOrd], name="useCov"),
    mxAlgebra(M[1,minOrd:maxOrd], name="useM"),
    mxFitFunctionML(),mxExpectationNormal(covariance="useCov", means="useM", dimnames = ordinalNames, threshnames = ordinalNames[1:nOrdinalVariables], thresholds="thresholdMatrix"),
    mxData(observed=jointData, type='raw')
)
checkPermutations(ordModel)

# Model spec for joint independent factors case
independentModel <- mxModel("independentComboModel",
    mxMatrix("Full", nVariables, nFactors, values=startLoads, free=startFree, lbound=-.99, ubound=.99, name="L"),
    mxMatrix("Unit", nVariables, 1, name="vectorofOnes"),
    mxMatrix("Zero", 1, nVariables, name="M"),
    mxAlgebra(vectorofOnes - (diag2vec(L %*% t(L))) , name="E"),
    mxAlgebra(L %*% t(L) + vec2diag(E), name="impliedCovs"),
    mxMatrix("Full", 
           name="thresholdDeviations", nrow=nThresholds, ncol=nOrdinalVariables,
           values=.2,
           free = useOptimizer, 
           lbound = rep( c(-Inf,rep(.01,(nThresholds-1))) , nOrdinalVariables),
           dimnames = list(c(), ordinalNames)),
   mxMatrix("Lower",nThresholds,nThresholds,values=1,free=F,name="unitLower"),
   mxAlgebra(unitLower %*% thresholdDeviations, name="thresholdMatrix"),
   mxFitFunctionML(),mxExpectationNormal(covariance="impliedCovs", means="M", dimnames = jointNames, threshnames = ordinalNames[1:nOrdinalVariables], thresholds="thresholdMatrix"),
   mxData(observed=jointData, type='raw')
)
checkPermutations(independentModel)

# Second Simulation: Continuous and Ordinal Variables no longer uncorrelated.
# We no longer have separability, so we'll run all-ordinal and all-continuous
# variants of the same model.
# For this, the joint results should fall between the ordinal and the continuous.

# Step 1: Set parameters for second simulation
loadings <- matrix(.4, nrow=nVariables,ncol=nFactors)  # Now, all correlated
startLoads <- loadings
residuals <- 1 - diag(loadings %*% t(loadings))
sigma <- loadings %*% t(loadings) + vec2diag(residuals)

mu <- matrix(0,nrow=nVariables,ncol=1)

# Step 2: simulate multivariate normal data for second simulation
set.seed(1234)
continuousCrossData <- data.frame(matrix(mvtnorm::rmvnorm(n=nSubjects,mu,sigma), nrow=nSubjects, ncol=nVariables))

# Step 3: chop continuous variables into ordinal data 
# with nThresholds+1 approximately equal categories, based on 1st variable
quants<-quantile(continuousCrossData[,1],  probs = c((1:nThresholds)/(nThresholds+1)))
ordinalCrossData<-matrix(0,nrow=nSubjects,ncol=nVariables)
for(i in 1:nVariables)
{
ordinalCrossData[,i] <- cut(as.vector(continuousCrossData[,i]),c(-Inf,quants,Inf))
}

# Step 4: make the ordinal variables into R factors and make a joint data frame with both variables in it, innit?
ordinalCrossData <- mxFactor(as.data.frame(ordinalCrossData),levels=c(1:(nThresholds+1)))
jointCrossData <- data.frame(continuousCrossData[,(nOrdinalVariables+1):nVariables],
	ordinalCrossData[,1:nOrdinalVariables])

names(jointCrossData) <- jointNames
names(continuousCrossData) <- jointNames
names(ordinalCrossData) <- jointNames

# Step 5: Set up actual models, simulation 2

# Ordinal-only model

thresholdModel <- mxModel("thresholdModel",
    mxMatrix("Full", nVariables, nFactors, values=startLoads, free=startFree, lbound=-.99, ubound=.99, name="L"),
    mxMatrix("Unit", nVariables, 1, name="vectorofOnes"),
    mxMatrix("Zero", 1, nVariables, name="M"),
    mxAlgebra(vectorofOnes - (diag2vec(L %*% t(L))) , name="E"),
    mxAlgebra(L %*% t(L) + vec2diag(E), name="impliedCovs"),
    mxMatrix("Full",
           name="thresholdDeviations", nrow=nThresholds, ncol=nVariables,
           values=.2,
           free = useOptimizer, 
           lbound = rep( c(-Inf,rep(.01,(nThresholds-1))) , nVariables),
           dimnames = list(c(), jointNames)),
    mxMatrix("Lower",nThresholds,nThresholds,values=1,free=F,name="unitLower"),
    mxAlgebra(unitLower %*% thresholdDeviations, name="thresholdMatrix"),
    mxFitFunctionML(),mxExpectationNormal(covariance="impliedCovs", means="M", dimnames = jointNames, thresholds="thresholdMatrix"),
           mxData(observed=ordinalCrossData, type='raw')
)

continuousModel <- mxModel("continuousModel",
    mxMatrix("Full", nVariables, nFactors, values=startLoads, free=startFree, lbound=-.99, ubound=.99, name="L"),
    mxMatrix("Unit", nVariables, 1, name="vectorofOnes"),
    mxMatrix("Zero", 1, nVariables, name="M"),
    mxAlgebra(vectorofOnes - (diag2vec(L %*% t(L))) , name="E"),
    mxAlgebra(L %*% t(L) + vec2diag(E), name="impliedCovs"),
    mxMatrix("Full", 
           name="thresholdDeviations", nrow=nThresholds, ncol=nVariables,
           values=.2,
           free = useOptimizer, 
           lbound = rep( c(-Inf,rep(.01,(nThresholds-1))) , nVariables),
           dimnames = list(c(), jointNames)),
    mxMatrix("Lower",nThresholds,nThresholds,values=1,free=F,name="unitLower"),
    mxAlgebra(unitLower %*% thresholdDeviations, name="thresholdMatrix"),
    mxFitFunctionML(),mxExpectationNormal(covariance="impliedCovs", means="M", dimnames = jointNames),
    mxData(observed=continuousCrossData, type='raw')
)

jointModel <- mxModel("jointModel",
    mxMatrix("Full", nVariables, nFactors, values=startLoads, free=startFree, lbound=-.99, ubound=.99, name="L"),
    mxMatrix("Unit", nVariables, 1, name="vectorofOnes"),
    mxMatrix("Zero", 1, nVariables, name="M"),
    mxAlgebra(vectorofOnes - (diag2vec(L %*% t(L))) , name="E"),
    mxAlgebra(L %*% t(L) + vec2diag(E), name="impliedCovs"),
    mxMatrix("Full", 
           name="thresholdDeviations", nrow=nThresholds, ncol=nOrdinalVariables,
           values=.2,
           free = useOptimizer, 
           lbound = rep( c(-Inf,rep(.01,(nThresholds-1))) , nOrdinalVariables),
           dimnames = list(c(), ordinalNames)),
   mxMatrix("Lower",nThresholds,nThresholds,values=1,free=F,name="unitLower"),
   mxAlgebra(unitLower %*% thresholdDeviations, name="thresholdMatrix"),
           mxFitFunctionML(),mxExpectationNormal(covariance="impliedCovs", means="M", dimnames = jointNames, threshnames = ordinalNames[1:nOrdinalVariables], thresholds="thresholdMatrix"),
           mxData(observed=jointCrossData, type='raw')
)

# Step 11: Set model options for speed

thresholdModel  <- mxOption(thresholdModel,  "Calculate Hessian", "No")
continuousModel <- mxOption(continuousModel, "Calculate Hessian", "No")
jointModel      <- mxOption(jointModel,      "Calculate Hessian", "No")
independentModel<- mxOption(independentModel,"Calculate Hessian", "No")
ordModel        <- mxOption(ordModel,        "Calculate Hessian", "No")
contModel       <- mxOption(contModel,       "Calculate Hessian", "No")
thresholdModel  <- mxOption(thresholdModel,  "Standard Errors", "No")
continuousModel <- mxOption(continuousModel, "Standard Errors", "No")
jointModel      <- mxOption(jointModel,      "Standard Errors", "No")
independentModel<- mxOption(independentModel,"Standard Errors", "No")
ordModel        <- mxOption(ordModel,        "Standard Errors", "No")
contModel       <- mxOption(contModel,       "Standard Errors", "No")

# Step 12: Run 'em.
thresholdModelRun <- mxRun(thresholdModel, suppressWarnings=TRUE)
continuousModelRun <- mxRun(continuousModel, suppressWarnings=TRUE)
jointModelRun <- mxRun(jointModel, suppressWarnings=TRUE)

ordinal    <- summary(thresholdModelRun)$Minus2LogLikelihood
continuous <- summary(continuousModelRun)$Minus2LogLikelihood
joint      <- summary(jointModelRun)$Minus2LogLikelihood

ord <- summary(mxRun(ordModel, suppressWarnings=TRUE))$Minus2LogLikelihood
cont <- summary(mxRun(contModel, suppressWarnings=TRUE))$Minus2LogLikelihood
indepJoint <- summary(mxRun(independentModel, suppressWarnings=TRUE))$Minus2LogLikelihood

# Second simulation set: 2 correlated continuous, one NA ordinal

useOptimizer <- TRUE

nSubjects <- 200
nThresh <- 1

cov1 <- matrix( c(2, .7, 0,
                 .7,  3, 0,
                  0,  0, 1), 3, 3)
start1 <- c(2, .7, 0, 3, 0, 1)
mu1 <-  c(15, -6, 0)
cNames1 <- c("C1", "C2", "Ona")
oNames1 <- c("O1", "O2", "Cna")

nVars <- dim(cov1)[1]
# Continuous and Ordinal vs joint: Continuous and ordinal data with an NA column of the other type
#   Comparison condition is to an NA column of the same type--with missingness, it's
#   all the same to the math.
set.seed(1234)
continuousData1 <- data.frame(matrix(mvtnorm::rmvnorm(n=nSubjects,mu1,cov1), nrow=nSubjects, ncol=nVars))
continuousData1[,nVars] <- as.numeric(NA)
names(continuousData1) <- cNames1



contModel1 <- mxModel("contModel1", 
    mxData(continuousData1, type="raw"),
    mxMatrix("Symm", nVars, nVars, values=start1, free=useOptimizer, name="Cov"),
    mxMatrix("Full", 1, nVars, values=rep(0, nVars), free=useOptimizer, name="Mean"),
    mxMatrix("Full", 1, 1, values=c(0), free=FALSE, name="Thresh"),
    mxFitFunctionML(),mxExpectationNormal(covariance="Cov", means="Mean", dimnames=cNames1)
)

contModel1A <- mxModel("contModel1A", 
    mxData(continuousData1, type="raw"),
    mxMatrix("Symm", nVars, nVars, values=start1, free=useOptimizer, name="Cov"),
    mxMatrix("Full", 1, nVars, values=rep(0, nVars), free=useOptimizer, name="Mean"),
    mxFitFunctionML(),mxExpectationNormal(covariance="Cov", means="Mean", dimnames=cNames1)
)

contModel1 <- mxOption(contModel1, "Function precision", 1e-9)
contModel1A <- mxOption(contModel1A, "Function precision", 1e-9)
contModel1 <- mxOption(contModel1, "Calculate Hessian", "No")
contModel1A <- mxOption(contModel1A, "Calculate Hessian", "No")
contModel1 <- mxOption(contModel1, "Standard Errors", "No")
contModel1A <- mxOption(contModel1A, "Standard Errors", "No")

contFit1 <- mxRun(contModel1, suppressWarnings=TRUE)
contFit1A <- mxRun(contModel1A, suppressWarnings=TRUE)
contSum1 <- summary(contFit1)
contSum1A <- summary(contFit1A)

ordinalData1 <- data.frame(matrix(mvtnorm::rmvnorm(n=nSubjects,mu1,cov1), nrow=nSubjects, ncol=nVars))
quants <- quantile(ordinalData1[,1],  probs = c((1:nThresh)/(nThresh+1)))
for(i in 1:nVars) {
   ordinalData1[,i] <- cut(as.vector(ordinalData1[,i]),c(-Inf,quants,Inf), labels=c(0:nThresh))
}
ordinalData1 <- mxFactor(ordinalData1, levels=c(0:nThresh))
ordinalData1[,nVars] <- mxFactor(NA, levels=1:2)
names(ordinalData1) <- oNames1
str(ordinalData1)

ordModel1 <- mxModel("ordModel1", 
    mxData(ordinalData1, type="raw"),
    mxMatrix("Symm", nVars, nVars, values=start1, free=useOptimizer, name="Cov"),
    mxMatrix("Full", 1, nVars, values=rep(0, nVars), free=useOptimizer, name="Mean"),
    mxMatrix("Full", nThresh, nVars, values=seq(-1, 1, length.out=nThresh), free=FALSE, name="Thresh"),
    mxFitFunctionML(),mxExpectationNormal(covariance="Cov", means="Mean", dimnames=oNames1, thresholds="Thresh", threshnames=oNames1)
)
# ordinalData1[,nVars] <- as.numeric(NA)
ordModel1A <- mxModel("ordModel1A", 
    mxData(ordinalData1, type="raw"),
    mxMatrix("Symm", nVars, nVars, values=start1, free=useOptimizer, name="Cov"),
    mxMatrix("Full", 1, nVars, values=rep(0, nVars), free=useOptimizer, name="Mean"),
    mxMatrix("Full", nThresh, nVars-1, values=seq(-1, 1, length.out=nThresh), free=FALSE, name="Thresh"),
    mxFitFunctionML(),mxExpectationNormal(covariance="Cov", means="Mean", dimnames = oNames1, thresholds="Thresh", threshnames=oNames1[1:(nVars-1)])
)

ordModel1 <- mxOption(ordModel1, "Function precision", 1e-9)
ordModel1A <- mxOption(ordModel1A, "Function precision", 1e-9)
ordModel1 <- mxOption(ordModel1, "Calculate Hessian", "No")
ordModel1A <- mxOption(ordModel1A, "Calculate Hessian", "No")
ordModel1 <- mxOption(ordModel1, "Standard Errors", "No")
ordModel1A <- mxOption(ordModel1A, "Standard Errors", "No")

ordFit1 <- mxRun(ordModel1, suppressWarnings=TRUE)
ordFit1A <- mxRun(ordModel1A, suppressWarnings=TRUE)
ordSum1 <- summary(ordFit1)
ordSum1A <- summary(ordFit1A)


# Second simulation: Same conditions as above, combined into a single model of two uncorrelated systems.
#   estimates should be identical.

nOrd <- 2
nCont <- 2
nVars <- nOrd + nCont
cov2 <- matrix(c(3, .7, 0, 0,
                 .7, 3, 0, 0,
                 0, 0, 1, -.7,
                 0, 0, -.7, 5), nVars,nVars)
mu2 <- c(5, 10, 0, 0)
startCont2 <- c(2, .7, 3)
startOrd2 <- c(1, .5, 1)
startAll2 <- c(1, .5, 0, 0, 1, 0, 0, 1, .5, 1)
startMeans2 <-c(0, 0, 0, 0)
allContinuousData2 <- data.frame(matrix(mvtnorm::rmvnorm(n=nSubjects,mu2,cov2), nrow=nSubjects, ncol=nVars))
continuousData2 <- allContinuousData2[,1:nCont]
ordinalData2 <- allContinuousData2[,(1:nOrd) + nCont]
quants <- quantile(ordinalData2[,1],  probs = c((1:nThresh)/(nThresh+1)))
for(i in 1:nOrd) {
   ordinalData2[,i] <- mxFactor(cut(as.vector(ordinalData2[,i]),c(-Inf,quants,Inf), labels=c(0:nThresh)), levels=0:nThresh)
}

cNames2 <- paste("C", 1:nCont, sep="")
oNames2 <- paste("O", 1:nOrd, sep="")
allNames2 <- c(cNames2, oNames2)
allData2 <- data.frame(continuousData2, ordinalData2)
names(continuousData2) <- cNames2
names(ordinalData2) <- oNames2
names(allData2) <- allNames2
contBound2 <- matrix(NA, nCont, nCont)
diag(contBound2) <- .01
ordBound2 <- matrix(NA, nOrd, nOrd)
diag(ordBound2) <- .01
allBound2 <- cbind(rbind(contBound2, matrix(NA, nOrd, nCont)), rbind(matrix(NA, nCont, nOrd), ordBound2))


contModel2 <- mxModel("contModel2", 
    mxData(continuousData2, type="raw"),
    mxMatrix("Symm", nCont, nCont, values=startCont2, free=useOptimizer, lbound=contBound2, name="Cov"),
    mxMatrix("Full", 1, nCont, values=startMeans2[1:nCont], free=useOptimizer, name="Mean"),
    mxFitFunctionML(),
    mxExpectationNormal(covariance="Cov", means="Mean", dimnames = cNames2)
    )

ordModel2 <-     mxModel("ordModel2", 
    mxData(ordinalData2, type="raw"),
    mxMatrix("Symm", nOrd, nOrd, values=startOrd2, free=useOptimizer, lbound=ordBound2, name="Cov"),
    mxMatrix("Full", 1, nCont, values=startMeans2[1:nOrd + nCont], free=useOptimizer, name="Mean"),
    mxMatrix("Full", nThresh, nOrd, values=seq(-1, 1, length.out=nThresh), free=FALSE, name="Thresh"),
    mxFitFunctionML(),
    mxExpectationNormal(covariance="Cov", means="Mean", dimnames = oNames2, thresholds="Thresh")
    )

allModel2 <- mxModel("jointModel2", 
    mxData(allData2, type="raw"),
    mxMatrix("Symm", nVars, nVars, values=startAll2, free=as.logical(useOptimizer * startAll2), lbound=allBound2, name="Cov"),
    mxMatrix("Full", 1, nVars, values=startMeans2[1:nOrd + nCont], free=useOptimizer, name="Mean"),
    mxMatrix("Full", nThresh, nOrd, values=seq(-1, 1, length.out=nThresh), free=FALSE, name="Thresh"),
    mxFitFunctionML(),
    mxExpectationNormal(covariance="Cov", means="Mean", dimnames = allNames2, thresholds="Thresh", threshnames = oNames2)
    )
    
dubData2 <- rbind(allData2, allData2)

dubModel2 <- mxModel("jointModel2Double", 
    mxData(dubData2, type="raw"),
    mxMatrix("Symm", nVars, nVars, values=startAll2, free=as.logical(useOptimizer * startAll2), lbound=allBound2, name="Cov"),
    mxMatrix("Full", 1, nVars, values=startMeans2[1:nOrd + nCont], free=useOptimizer, name="Mean"),
    mxMatrix("Full", nThresh, nOrd, values=seq(-1, 1, length.out=nThresh), free=FALSE, name="Thresh"),
    mxFitFunctionML(),mxExpectationNormal(covariance="Cov", means="Mean", dimnames = allNames2, thresholds="Thresh", threshnames = oNames2)
    )


contModel2 <- mxOption(contModel2, "Function precision", 1e-9)
ordModel2 <- mxOption(ordModel2, "Function precision", 1e-9)
allModel2 <- mxOption(allModel2, "Function precision", 1e-9)
dubModel2 <- mxOption(dubModel2, "Function precision", 1e-9)

contModel2 <- mxOption(contModel2, "Calculate Hessian", "No")
ordModel2 <- mxOption(ordModel2, "Calculate Hessian", "No")
allModel2 <- mxOption(allModel2, "Calculate Hessian", "No")
dubModel2 <- mxOption(dubModel2, "Calculate Hessian", "No")

contModel2 <- mxOption(contModel2, "Standard Errors", "No")
ordModel2 <- mxOption(ordModel2, "Standard Errors", "No")
allModel2 <- mxOption(allModel2, "Standard Errors", "No")
dubModel2 <- mxOption(dubModel2, "Standard Errors", "No")

contFit2 <- mxRun(contModel2, suppressWarnings=TRUE)
ordFit2  <- mxRun(ordModel2, suppressWarnings=TRUE)
allFit2  <- mxRun(allModel2, suppressWarnings=TRUE)
dubFit2  <- mxRun(dubModel2, suppressWarnings=TRUE)

contSum2 <- summary(contFit2)
ordSum2  <- summary(ordFit2)
allSum2  <- summary(allFit2)
dubSum2  <- summary(dubFit2)

# Second simulation: The same conditions from above, with correlations freed between.
#   Need correct values to make a solid check
# 
# nOrd <- 2
# nCont <- 2
# nVars <- nOrd + nCont
# cov3 <- matrix(c(3, .7, .8, -.2,
#                  .7, 3, .6, -.4,
#                  .8, .6, 1, -.7,
#                  -.2, -.4, -.7, 5), nVars,nVars)
# mu3 <- c(5, 10, 0, 0)
# startAll3 <- c(1, .5, .5, .5, 1, .5, .5, 1, .5, 1)
# startMeans3 <-c(0, 0, 0, 0)
# allContinuousData3 <- data.frame(matrix(mvtnorm::rmvnorm(n=nSubjects,mu3,cov3), nrow=nSubjects, ncol=nVars))
# continuousData3 <- allContinuousData3[,1:nCont]
# ordinalData3 <- allContinuousData3[,(1:nOrd) + nCont]
# quants <- quantile(ordinalData3[,1],  probs = c((1:nThresh)/(nThresh+1)))
# for(i in 1:nOrd) {
#    ordinalData3[,i] <- mxFactor(cut(as.vector(ordinalData3[,i]),c(-Inf,quants,Inf), labels=c(0:nThresh)), levels=0:nThresh)
# }
# 
# cNames3 <- paste("C", 1:nCont, sep="")
# oNames3 <- paste("O", 1:nOrd, sep="")
# allNames3 <- c(cNames3, oNames3)
# allData3 <- data.frame(continuousData3, ordinalData3)
# names(allData3) <- allNames3
# contBound3 <- matrix(NA, nCont, nCont)
# diag(contBound3) <- .01
# ordBound3 <- matrix(NA, nOrd, nOrd)
# diag(ordBound3) <- .01
# allBound3 <- cbind(rbind(contBound3, matrix(NA, nOrd, nCont)), rbind(matrix(NA, nCont, nOrd), ordBound3))
# 
# allModel3 <- mxModel("jointModel3", 
#     mxData(allData3, type="raw"),
#     mxMatrix("Symm", nVars, nVars, values=startAll3, free=as.logical(useOptimizer * startAll3), lbound=allBound3, name="Cov"),
#     mxMatrix("Full", 1, nVars, values=startMeans3[1:nOrd + nCont], free=useOptimizer, name="Mean"),
#     mxMatrix("Full", nThresh, nOrd, values=seq(-1, 1, length.out=nThresh), free=FALSE, name="Thresh"),
#     mxFitFunctionML(),mxExpectationNormal(covariance="Cov", means="Mean", dimnames = allNames3, thresholds="Thresh", threshnames = oNames3)
#     )
# 
# allModel3 <- mxOption(allModel3, "Function precision", 1e-9)
# allModel3 <- mxOption(allModel3, "Calculate Hessian", "No")
# 
# allFit3 <- mxRun(allModel3, suppressWarnings=TRUE)
# allSum3 <- summary(allFit3)

omxCheckCloseEnough(contSum1$Minus2LogLikelihood, contSum1A$Minus2LogLikelihood, .0000001)
omxCheckCloseEnough(ordSum1$Minus2LogLikelihood, ordSum1A$Minus2LogLikelihood, .0000001)
omxCheckCloseEnough(contSum2$Minus2LogLikelihood + ordSum2$Minus2LogLikelihood, allSum2$Minus2LogLikelihood, .01)
omxCheckCloseEnough(2* allSum2$Minus2LogLikelihood, dubSum2$Minus2LogLikelihood, .01)

# Step 13: Tests

# Isbetween will be negative if the signs don't match; that is, if joint is not between
isbetween <- sign(ordinal - joint) * sign(joint - continuous)

omxCheckEquals(isbetween, 1)
omxCheckCloseEnough(indepJoint, ord+cont, .000001)




