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
# Program: OneFactorOrdinal01_MatrixRaw.R  
# Author: Michael Neale
# Date: 2010.08.14 
#
# ModelType: Factor
# DataType: Ordinal
# Field: None
#
# Purpose: 
#      One Factor model to estimate factor loadings, residual variances and means
#      Matrix style model input - Raw data input - Ordinal data
#
# RevisionHistory:
#      Hermine Maes -- 2010.09.07 updated & reformatted & respecified
#      Ross Gore -- 2011.06.06	added Model, Data & Field metadata
# -----------------------------------------------------------------------------

require(OpenMx)
require(MASS)
# Load Libraries
# -----------------------------------------------------------------------------

nVariables<-5
nFactors<-1
nThresholds<-3
nSubjects<-500
isIdentified<-function(nVariables,nFactors) as.logical(1+sign((nVariables*(nVariables-1)/2) - nVariables*nFactors + nFactors*(nFactors-1)/2))

isIdentified(nVariables,nFactors) # if this function returns FALSE then model is not identified, otherwise it is.
# Set up simulation parameters: 
# nVariables>=3, nThresholds>=1, 
# nSubjects>=nVariables*nThresholds 
# and model should be identified
# -------------------------------------

loadings <- matrix(.7,nrow=nVariables,ncol=nFactors)
residuals <- 1 - (loadings * loadings)
sigma <- loadings %*% t(loadings) + vec2diag(residuals)
mu <- matrix(0,nrow=nVariables,ncol=1)

set.seed(1234)
continuousData <- mvrnorm(n=nSubjects,mu,sigma)
# Simulate multivariate normal data
# -------------------------------------

quants<-quantile(continuousData[,1],  probs = c((1:nThresholds)/(nThresholds+1)))
ordinalData<-matrix(0,nrow=nSubjects,ncol=nVariables)
for(i in 1:nVariables)
{
ordinalData[,i] <- cut(as.vector(continuousData[,i]),c(-Inf,quants,Inf))
}
# Chop continuous variables into 
# ordinal data with nThresholds+1 
# approximately equal categories, 
# based on 1st variable
# -------------------------------------

ordinalData <- mxFactor(as.data.frame(ordinalData),levels=c(1:(nThresholds+1)))
# Make the ordinal variables into 
# R factors
# -------------------------------------

fruitynames<-paste("banana",1:nVariables,sep="")
names(ordinalData)<-fruitynames
# Name the variables
# -------------------------------------

# Simulate Data
# -----------------------------------------------------------------------------

oneFactorThresholdModel01 <- mxModel("oneFactorThresholdModel01",
    mxMatrix(
        type="Full", 
        nrow=nVariables, 
        ncol=nFactors, 
        free=TRUE, 
        values=0.2, 
        lbound=-.99, 
        ubound=2, 
        name="facLoadings"
    ),
    mxMatrix(
        type="Diag", 
        nrow=nVariables, 
        ncol=nVariables,
        free=TRUE,
        values=0.9,
        name="resVariances"
    ),
    mxAlgebra(
        expression=facLoadings %*% t(facLoadings) + resVariances, 
        name="expCovariances"
    ),
    mxMatrix(
        type="Full", 
        nrow=1, 
        ncol=nVariables,
        free=T,
        name="expMeans"
    ),
    mxMatrix(
        type="Full", 
        nrow=nThresholds, 
        ncol=nVariables,
        free=rep( c(F,F,rep(T,(nThresholds-2))), nVariables), 
        values=rep( c(0,1,rep(.2,(nThresholds-2))), nVariables),
        lbound=rep( c(-Inf,rep(.01,(nThresholds-1))), nVariables),
        dimnames=list(c(), fruitynames),
        name="thresholdDeviations"
    ),
    mxMatrix(
        type="Lower",
        nrow=nThresholds,
        ncol=nThresholds,
        free=FALSE,
        values=1,
        name="unitLower"
    ),
    mxAlgebra(
        expression=unitLower %*% thresholdDeviations, 
        name="expThresholds"
    ),
        mxMatrix(
    	 type="Unit",
    	 nrow=nThresholds,
    	 ncol=1,
    	 name="columnofOnes"
    ),
    mxAlgebra(
    	 expression=expMeans %x% columnofOnes,
    	 name="meansMatrix"
    ),
    mxAlgebra(
    	 expression=sqrt(t(diag2vec(expCovariances))) %x% columnofOnes,
    	 name="variancesMatrix"
    ),
    mxAlgebra(
    	 expression=(expThresholds - meansMatrix) / variancesMatrix,
    	 name="thresholdMatrix"
    ),
    mxMatrix( 
    	type="Iden", 
    	nrow=nVariables, 
    	ncol=nVariables, 
    	name="Identity"
    ),
    mxAlgebra(
    	 expression=solve(sqrt(Identity * expCovariances)) %*% facLoadings,
    	 name="facLoadingsMatrix"
    ),
    mxData(
        observed=ordinalData, 
        type='raw'
    ),
    mxFIMLObjective(
        covariance="expCovariances", 
        means="expMeans", 
        dimnames=fruitynames, 
        thresholds="expThresholds"
    )
)
# Create Factor Model with Raw Ordinal Data and Matrices Input
# -----------------------------------------------------------------------------

oneFactorThresholdFit01 <- mxRun(oneFactorThresholdModel01, suppressWarnings=TRUE)
# Fit the model with mxRun
# -----------------------------------------------------------------------------

summary(oneFactorThresholdFit01)
# Print a summary of the results
# -----------------------------------------------------------------------------