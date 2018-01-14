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

# -----------------------------------------------------------------------------
# Program: OneFactorOrdinal_MatrixRaw.R  
# Author: Michael Neale
# Date: 2010.08.14 
#
# ModelType: Factor
# DataType: Ordinal
# Field: None
#
# Purpose:
#      One Factor model to estimate factor loadings, 
#      residual variances and means
#      Matrix style model input - Raw data input - Ordinal data
#
# RevisionHistory:
#      Hermine Maes -- 2010.09.07 updated & reformatted
#      Ross Gore -- 2011.06.06	added Model, Data & Field metadata
# -----------------------------------------------------------------------------

require(OpenMx)

# Load libraries
# -----------------------------------------------------------------------------

nVariables<-5
nFactors<-1
nThresholds<-3
nSubjects<-500
isIdentified<-function(nVariables,nFactors) as.logical(1+sign((nVariables*(nVariables-1)/2) -  nVariables*nFactors + nFactors*(nFactors-1)/2))

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
continuousData <- mvtnorm::rmvnorm(n=nSubjects,mu,sigma)
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

oneFactorThresholdModel <- mxModel("oneFactorThresholdModel",
    mxMatrix(
        type="Full", 
        nrow=nVariables, 
        ncol=nFactors, 
        free=TRUE, 
        values=0.2, 
        lbound=-.99, 
        ubound=.99, 
        name="facLoadings"
    ),
    mxMatrix(
        type="Unit", 
        nrow=nVariables, 
        ncol=1, 
        name="vectorofOnes"
    ),
    mxAlgebra(
        expression=vectorofOnes - (diag2vec(facLoadings %*% t(facLoadings))) , 
        name="resVariances"
    ),
    mxAlgebra(
        expression=facLoadings %*% t(facLoadings) + vec2diag(resVariances), 
        name="expCovariances"
    ),
    mxMatrix(
        type="Zero", 
        nrow=1, 
        ncol=nVariables, 
        name="expMeans"
    ),
    mxMatrix(
        type="Full", 
        nrow=nThresholds, 
        ncol=nVariables,
        free=TRUE, 
        values=.2,
        lbound=rep( c(-Inf,rep(.01,(nThresholds-1))) , nVariables),
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
    mxData(
        observed=ordinalData, 
        type='raw'
    ),
    mxFitFunctionML(),mxExpectationNormal(
        covariance="expCovariances", 
        means="expMeans", 
        dimnames=fruitynames, 
        thresholds="expThresholds"
    )
)
# Create Factor Model with Raw Ordinal Data and Matrices Input
# -----------------------------------------------------------------------------

oneFactorThresholdFit <- mxRun(oneFactorThresholdModel, suppressWarnings=TRUE)
# Fit the model with mxRun
# -----------------------------------------------------------------------------

summary(oneFactorThresholdFit)
# Print a summary of the results
# -----------------------------------------------------------------------------

omxCheckCloseEnough(oneFactorThresholdFit$output$fit, 6281.996, .1)
