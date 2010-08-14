#
#  OpenMx Ordinal Data Example
#  Revision history:
#		Michael Neale 14 Aug 2010
#

# Step 1: load libraries
require(OpenMx)
require(MASS)

# Step 2: set up simulation parameters
nVariables<-8
nFactors<-1
nThresholds<-3
nSubjects<-500

loadings <- matrix(.7,nrow=nVariables,ncol=nFactors)
residuals <- 1 - (loadings * loadings)
sigma <- loadings %*% t(loadings) + vec2diag(residuals)
mu <- matrix(0,nrow=nVariables,ncol=1)
# Step 3: simulate multivariate normal data
continuousData <- mvrnorm(n=nSubjects,mu,sigma)

# Step 4: chop continuous variables into ordinal data
ordinalData<-matrix(0,nrow=nSubjects,ncol=nVariables)
for(i in 1:nVariables)
{
ordinalData[,i] <- cut(as.vector(continuousData[,i]),c(-Inf,-1,0,1,Inf))
}

# Step 5: make the ordinal variables into R factors
ordinalData <- mxFactor(as.data.frame(ordinalData),levels=c(1,2,3,4))

# Step 6: name the variables
fruitynames<-paste("banana",1:nVariables,sep="")
names(ordinalData)<-fruitynames


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
            mxFIMLObjective(covariance="impliedCovs", means="M", dimnames = fruitynames, thresholds="thresholdMatrix"),
            mxData(observed=ordinalData, type='raw')
)

summary(thresholdModelrun <- mxRun(thresholdModel))
#
#  Now try it with the 0/1 parameterization and standardization
#
thresholdModel2 <- mxModel("thresholdModel",
	mxMatrix("Full", nVariables, nFactors, values=0.2, free=TRUE, name="L"),
	mxMatrix("Unit", nThresholds, 1, name="U"),
	mxMatrix("Full", 1, nVariables, name="M", free=TRUE),
	mxMatrix("Diag", nVariables, nVariables, name="E", free=TRUE),
	mxAlgebra(L %*% t(L) + E, name="impliedCovs"),
	mxMatrix("Full", 
            name="thresholdDeviations",
            values=rbind(rep(0,nVariables),rep(1,nVariables),rep(.1,((nThresholds-2)*nVariables))),
            free = rbind(rep(FALSE, nVariables),rep(FALSE, nVariables),rep(TRUE, ((nThresholds-2)*nVariables))),
            lbound = rep( c(-Inf,rep(.001,(nThresholds-1))) , nVariables),
            dimnames = list(c(), fruitynames)),
    mxMatrix("Lower",nThresholds,nThresholds,values=1,free=F,name="unitLower"),
    mxAlgebra(unitLower %*% thresholdDeviations, name="thresholdMatrix"),
            mxFIMLObjective(covariance="impliedCovs", means="M", dimnames = fruitynames, thresholds="thresholdMatrix"),
            mxData(observed=ordinalData, type='raw')
)

summary(thresholdModel2run <- mxRun(thresholdModel2))

# # Now we standardize the 0/1 fixed threshold estimates # using a nthresholds x 1 Unit vector U to reproduce rows etc. #
standardizedThresholds <- mxEval((thresholdMatrix - ((U %x% M[,1:3])) ) / (U %x% t(sqrt(diag2vec(impliedCovs)))) , thresholdModel2run)
# Flushed with success at that, we now go after Standardized A and S matrices (loadings & initial covariances)
D <- mxEval(sqrt(diag2vec(impliedCovs)), thresholdModel2run)
standardizedLoadings <- mxEval( L /  D, thresholdModel2run)
# These are reasonably close to those from thresholdModel:
thresholdModelrun$L@values - standardizedLoadings
# as are the standardized thresholds
standardizedThresholds - thresholdModelrun$thresholdMatrix@result