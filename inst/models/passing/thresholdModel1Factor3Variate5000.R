# History
# 2017-04-14 04:28PM TBATES update for mxExpectationNormal

# Step 1: load libraries
require(OpenMx)
require(MASS)
#
# Step 2: set up simulation parameters 
# Note: nVariables>=3, nThresholds>=1, nSubjects>=nVariables*nThresholds (maybe more)
# and model should be identified
#

nVariables  <- 3
nFactors    <- 1
nThresholds <- 3
nSubjects   <- 5000
isIdentified<-function(nVariables,nFactors) as.logical(1+sign((nVariables*(nVariables-1)/2) -  nVariables*nFactors + nFactors*(nFactors-1)/2))
# if this function returns FALSE then model is not identified, otherwise it is.
isIdentified(nVariables,nFactors)

loadings  <- matrix(.7,nrow=nVariables,ncol=nFactors)
residuals <- 1 - (loadings * loadings)
sigma     <- loadings %*% t(loadings) + vec2diag(residuals)
mu        <- matrix(0,nrow=nVariables,ncol=1)
# Step 3: simulate multivariate normal data
set.seed(1234)
continuousData <- mvrnorm(n=nSubjects,mu,sigma)

# Step 4: chop continuous variables into ordinal data 
# with nThresholds+1 approximately equal categories, based on 1st variable
quants      <- quantile(continuousData[,1], probs = c((1:nThresholds)/(nThresholds+1)))
ordinalData <- matrix(0,nrow=nSubjects,ncol=nVariables)
for(i in 1:nVariables){
	ordinalData[,i] <- cut(as.vector(continuousData[,i]),c(-Inf,quants,Inf))
}

# Step 5: make the ordinal variables into R factors
ordinalData <- mxFactor(as.data.frame(ordinalData),levels=c(1:(nThresholds+1)))

# Step 6: name the variables
fruitynames <- paste("banana",1:nVariables,sep="")
names(ordinalData) <- fruitynames


thresholdModel <- mxModel("thresholdModel",
	mxMatrix("Full", nVariables, nFactors, values=0.2, free=TRUE, lbound=-.99, ubound=.99, name="L"),
	mxMatrix("Unit", nVariables, 1, name="vectorofOnes"),
	mxMatrix("Zero", 1, nVariables, name="M"),
	mxAlgebra(vectorofOnes - (diag2vec(L %*% t(L))) , name="E"),
	mxAlgebra(L %*% t(L) + vec2diag(E), name="impliedCovs"),
	mxMatrix("Full", 
            name="thresholdDeviations", nrow=nThresholds, ncol=nVariables,
            values=.2,
            free = TRUE, 
            lbound = rep( c(-Inf,rep(.01,(nThresholds-1))) , nVariables),
            dimnames = list(c(), fruitynames)),
		mxMatrix("Lower",nThresholds,nThresholds,values=1,free=F,name="unitLower"),
		mxAlgebra(unitLower %*% thresholdDeviations, name="thresholdMatrix"),
		mxExpectationNormal(covariance="impliedCovs", means="M", dimnames = fruitynames, thresholds="thresholdMatrix"),
    mxFitFunctionML(),
		mxData(observed=ordinalData, type='raw')
)

thresholdModelrun <- mxRun(thresholdModel)
# now (OpenMx version: 2.7.9.58) giving code Red
# In model 'thresholdModel' Optimizer returned a non-zero status code 6. The model does not satisfy the first-order optimality conditions to the required accuracy, and no improved point for the merit function could be found during the final linesearch (Mx status RED)

summary(thresholdModelrun)

# NOTE! The reason that this ____USED TO BE____ a failing test is because Windows can
# find a better solution, 39011.83. Possibly there are some box
# constraints missing.
omxCheckCloseEnough(thresholdModelrun$output$fit, 39034.359, 1.5)

got <- omxGetParameters(thresholdModelrun)
names(got) <- NULL
est <- c(0.693, 0.73, 0.695, -0.674, 0.678, 0.672, -0.626,  0.639, 0.642, -0.649, 0.659, 0.653)
omxCheckCloseEnough(got, est, .02)
