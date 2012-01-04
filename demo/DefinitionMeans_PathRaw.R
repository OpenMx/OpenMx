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

# -----------------------------------------------------------------------------
# Program: DefinitionMeans_PathRaw.R  
# Author: Mike Neale
# Date: 2009.08.01 
#
# ModelType: Means
# DataType: Continuous
# Field: None
#
# Purpose: 
#      Definition Means model to estimate moderation effect 
#      of measured variable 
#      Path style model input - Raw data input
#
# RevisionHistory:
#      Hermine Maes -- 2009.10.08 updated & reformatted
#      Ross Gore -- 2011.06.15 added Model, Data & Field metadata
# -----------------------------------------------------------------------------

#This script is used to test the definition variable functionality in OpenMx.
#The definition variable in this example is dichotomous, and describes 
# two different groups.
#These two groups are measured on two variables, x and y.
#The group with a definition value of 1 has means of 1 and 2 for x and y.
#The group with a definition value of 0 has means af zero for x and y. 
#The definition variable is used to define a mean deviation of the group
# with definition value 1.

require(OpenMx)
library(MASS)
# Load Libraries
# -----------------------------------------------------------------------------

set.seed(200)
N=500
Sigma <- matrix(c(1,.5,.5,1),2,2)
group1<-mvrnorm(N, c(1,2), Sigma) # use mvrnorm from MASS package
group2<-mvrnorm(N, c(0,0), Sigma)

y<-rbind(group1,group2)           # Bind both groups together by rows
dimnames(y)[2]<-list(c("x","y")); # Add names
def    <-rep(c(1,0),each=N);      # Add a definition variable 2n 
								  # in length for group status
selVars<-c("x","y")               # Make a selection variables object
# Simulate data
# -----------------------------------------------------------------------------

defMeansModel <- mxModel("Definition Means Path Specification", 
	type="RAM",
	manifestVars=selVars,
	latentVars  ="DefDummy",
    mxPath(
    	from=c("x","y"), 
    	arrows=2, 
    	free= TRUE, 
    	values=1,  
    	labels=c("Varx","Vary")
    ),
	# variances
	# -------------------------------------
    mxPath(
    	from="x", 
    	to="y", 
    	arrows=2, 
    	free= TRUE, 
    	values=.1, 
    	labels=c("Covxy")
    ),
	# covariances
	# -------------------------------------        
    mxPath(
    	from="one", 
    	to=c("x","y"), 
    	arrows=1, 
    	free= TRUE, 
    	values=1, 
    	labels=c("meanx","meany")
    ),
 	# means 
	# ------------------------------------- 
    mxPath(
    	from="one", 
    	to="DefDummy", 
    	arrows=1, 
    	free= FALSE, 
    	values=1, 
    	labels="data.def"
    ),
    # definition value
	# ------------------------------------- 
    mxPath(
    	from="DefDummy", 
    	to=c("x","y"), 
    	arrows=1, 
    	free= TRUE, 
    	values=1, 
    	labels=c("beta_1","beta_2")
    ), 
	# beta weights
	# -------------------------------------
    mxData(
    	observed=data.frame(y,def), 
    	type="raw"
    )
)
#Define model
# -----------------------------------------------------------------------------

defMeansFit<-mxRun(defMeansModel)
# Run the model
# -----------------------------------------------------------------------------

defMeansFit@matrices
defMeansFit@algebras



# Remember to knock off 1 and 2 
# from group 1's data
# so as to estimate variance of 
# combined sample without the mean 
# correction. First we compute some 
# summary statistics from the data
# -------------------------------------
ObsCovs        <- cov(rbind(group1 - rep(c(1,2), each=N), group2))
ObsMeansGroup1 <- c(mean(group1[,1]), mean(group1[,2]))
ObsMeansGroup2 <- c(mean(group2[,1]), mean(group2[,2]))

# Second we extract the parameter 
# estimates and matrix algebra results 
# from the model.
# -------------------------------------
Sigma <- mxEval(S[1:2,1:2], defMeansFit)
Mu    <- mxEval(M[1:2], defMeansFit)
beta  <- mxEval(A[1:2,3], defMeansFit)

# Third, we check to see if things are
# more or less equal.
# -------------------------------------
omxCheckCloseEnough(ObsCovs,Sigma,.01)
omxCheckCloseEnough(ObsMeansGroup1,as.vector(Mu+beta),.001)
omxCheckCloseEnough(ObsMeansGroup2,as.vector(Mu),.001)
# Compare OpenMx estimates to summary statistics from raw data, 
# -----------------------------------------------------------------------------
