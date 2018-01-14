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
# Program: DefinitionMeans_MatrixRaw.R  
# Author: Mike Neale
# Date: 2009.08.01 
#
# ModelType: Means
# DataType: Continuous
# Field: None
#
# Purpose: 
#      Definition Means model to estimate moderation effect of measured variable 
#      Matrix style model input - Raw data input
#
# RevisionHistory:
#      Hermine Maes -- 2009.10.08 updated & reformatted
#      Ross Gore -- 2011.06.15 added Model, Data & Field
#      Hermine Maes -- 2014.11.02 piecewise specification
# -----------------------------------------------------------------------------

#This script is used to test the definition variable functionality in OpenMx
#The definition variable in this example is dichotomous, and describes two different groups
#These two groups are measured on two variables, x and y
#The group with a definition value of 1 has means of 1 and 2 for x and y
#The group with a definition value of 0 has means af zero for x and y 
#The definition variable is used to define a mean deviation of the group with definition value 1


require(OpenMx)
library(MASS)
# Load Library
# -----------------------------------------------------------------------------

set.seed(200)
N=500
Sigma          <- matrix(c(1,.5,.5,1),2,2)
group1         <- mvtnorm::rmvnorm(N, c(1,2), Sigma) # Use mvrnorm from MASS package
group2         <- mvtnorm::rmvnorm(N, c(0,0), Sigma)

xy             <- rbind(group1,group2)      # Bind groups together by rows
dimnames(xy)[2]<- list(c("x","y"))          # Add names
def            <- rep(c(1,0),each=N);       # Add def var [2n] for group status
selVars        <- c("x","y")                # Make selection variables object
# Simulate data
# -----------------------------------------------------------------------------

Sigma        <- mxMatrix( type="Symm", nrow=2, ncol=2, 
    	                  free=TRUE, values=c(1, 0, 1), name="Sigma" )
Mean         <- mxMatrix( type="Full", nrow=1, ncol=2, 
                          free=TRUE, name="Mean" )
beta         <- mxMatrix( type="Full", nrow=1, ncol=2, 
	                      free=TRUE, values=c(0,0), name="beta" )
dataDef      <- mxMatrix( type="Full", nrow=1, ncol=2, 
    	                  free=FALSE, labels=c("data.def"), name="def" )
Mu           <- mxAlgebra( expression=Mean + beta*def, name="Mu" )
dataRaw      <- mxData( observed=data.frame(xy,def), type="raw" )
exp          <- mxExpectationNormal( covariance="Sigma", means="Mu", 
                                     dimnames=selVars )
funML        <- mxFitFunctionML()
defMeansModel <- mxModel("Definition  Means Matrix Specification", 
                         dataRaw, Sigma, Mean, beta, dataDef, Mu, exp, funML)
# Define model
# -----------------------------------------------------------------------------


defMeansFit <- mxRun(defMeansModel)
# Run the model
# -----------------------------------------------------------------------------

defMeansFit$matrices
defMeansFit$algebras



ObsCovs <- cov(rbind(group1 - rep(c(1,2), each=N), group2))
ObsMeansGroup1 <- c(mean(group1[,1]), mean(group1[,2]))
ObsMeansGroup2 <- c(mean(group2[,1]), mean(group2[,2]))
# Remember to knock off 1 and 2 from group 1's data
# so as to estimate variance of combined 
# sample without the mean correction.
# First we compute some summary 
# statistics from the data
# -------------------------------------

Sigma <- mxEval(Sigma, defMeansFit)
Mean <- mxEval(Mean, defMeansFit)
beta <- mxEval(beta, defMeansFit)
# Second we extract the parameter 
# estimates and matrix algebra results 
# from the model
# -------------------------------------

omxCheckCloseEnough(ObsCovs, Sigma, .01)
omxCheckCloseEnough(ObsMeansGroup1, as.vector(Mean+beta), .001)
omxCheckCloseEnough(ObsMeansGroup2, as.vector(Mean), .001)
# Third, we check to see if things are 
# more or less equal
# -------------------------------------
# Compare OpenMx estimates to summary statistics from raw data,
# -----------------------------------------------------------------------------
