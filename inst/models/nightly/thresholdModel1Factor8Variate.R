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
#  OpenMx Ordinal Data Example
#  Revision history:
#		Michael Neale 14 Aug 2010
#

# Step 1: load libraries
require(OpenMx)

#
# Step 2: set up simulation parameters 
# Note: nVariables>=3, nThresholds>=1, nSubjects>=nVariables*nThresholds (maybe more)
# and model should be identified
#
nVariables<-8
nFactors<-1
nThresholds<-3
nSubjects<-500
isIdentified<-function(nVariables,nFactors) as.logical(1+sign((nVariables*(nVariables-1)/2) -  nVariables*nFactors + nFactors*(nFactors-1)/2))
# if this function returns FALSE then model is not identified, otherwise it is.
isIdentified(nVariables,nFactors)

loadings <- matrix(.7,nrow=nVariables,ncol=nFactors)
residuals <- 1 - (loadings * loadings)
sigma <- loadings %*% t(loadings) + vec2diag(residuals)
mu <- matrix(0,nrow=nVariables,ncol=1)
# Step 3: simulate multivariate normal data
set.seed(1234)
continuousData <- mvtnorm::rmvnorm(n=nSubjects,mu,sigma)

# Step 4: chop continuous variables into ordinal data 
# with nThresholds+1 approximately equal categories, based on 1st variable
quants<-quantile(continuousData[,1],  probs = c((1:nThresholds)/(nThresholds+1)))
ordinalData<-matrix(0,nrow=nSubjects,ncol=nVariables)
for(i in 1:nVariables)
{
ordinalData[,i] <- cut(as.vector(continuousData[,i]),c(-Inf,quants,Inf))
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
            mxFitFunctionML(),mxExpectationNormal(covariance="impliedCovs", means="M", dimnames = fruitynames, thresholds="thresholdMatrix"),
            mxData(observed=ordinalData, type='raw')
)

summary(thresholdModelrun <- mxTryHardOrdinal(thresholdModel, extraTries=1L))
omxCheckCloseEnough(thresholdModelrun$output$fit, 9764.721, .1)
