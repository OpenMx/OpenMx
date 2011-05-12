#
#  OpenMx Mixed Ordinal & Continuous Data Example
#  Revision history:
#		Michael Neale 7 Jan 2011
#       Revised Timothy Brick 10 Jan 2011

require(OpenMx)
require(MASS)

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
continuousData <- data.frame(matrix(mvrnorm(n=nSubjects,mu,sigma), nrow=nSubjects, ncol=nVariables))

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
    mxFIMLObjective(covariance="useCov", means="useM", dimnames = continuousNames),
    mxData(observed=jointData, type='raw')
)

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
    mxFIMLObjective(covariance="useCov", means="useM", dimnames = ordinalNames, threshnames = ordinalNames[1:nOrdinalVariables], thresholds="thresholdMatrix"),
    mxData(observed=jointData, type='raw')
)

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
   mxFIMLObjective(covariance="impliedCovs", means="M", dimnames = jointNames, threshnames = ordinalNames[1:nOrdinalVariables], thresholds="thresholdMatrix"),
   mxData(observed=jointData, type='raw')
)

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
continuousCrossData <- data.frame(matrix(mvrnorm(n=nSubjects,mu,sigma), nrow=nSubjects, ncol=nVariables))

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
    mxFIMLObjective(covariance="impliedCovs", means="M", dimnames = jointNames, thresholds="thresholdMatrix"),
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
    mxFIMLObjective(covariance="impliedCovs", means="M", dimnames = jointNames),
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
           mxFIMLObjective(covariance="impliedCovs", means="M", dimnames = jointNames, threshnames = ordinalNames[1:nOrdinalVariables], thresholds="thresholdMatrix"),
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
thresholdModelRun <- mxRun(thresholdModel)
continuousModelRun <- mxRun(continuousModel)
jointModelRun <- mxRun(jointModel)

ordinal    <- summary(thresholdModelRun)$Minus2LogLikelihood
continuous <- summary(continuousModelRun)$Minus2LogLikelihood
joint      <- summary(jointModelRun)$Minus2LogLikelihood

ord <- summary(mxRun(ordModel))$Minus2LogLikelihood
cont <- summary(mxRun(contModel) )$Minus2LogLikelihood
indepJoint <- summary(mxRun(independentModel))$Minus2LogLikelihood

# Step 13: Tests

# Isbetween will be negative if the signs don't match; that is, if joint is not between
isbetween <- sign(ordinal - joint) * sign(joint - continuous)

omxCheckEquals(isbetween, 1)
omxCheckCloseEnough(indepJoint, ord+cont, .000001)




