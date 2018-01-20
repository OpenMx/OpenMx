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

require(OpenMx)

#Latent class analysis example 
#------------------------------------------------------+
# Marginal ml                                          |
# data from http://www.people.vcu.edu/~nhenry/LSA50.htm|
# see                                                  |
# Latent Class Probs Item 1 Item 2 Item 3 Item 4       |
# .4243              .9240 .6276 .5704 .5125           |
# .5757              .4324 .1871 .1008 .0635           |
# Mx recovers:                                         |
# .4440    0.90736    0.61444    0.55828    0.50176    |
# .5560    0.42838    0.18218    .09394      .05627    |
#------------------------------------------------------+

# Data
cn <- c("Armyrun", "Favatt", "squaredeal", "welfare", "freq")
data <- suppressWarnings(try(read.table("models/passing/data/lazarsfeld.ord", col.names=cn), silent=TRUE))
if (is(data, "try-error")) data <- read.table("data/lazarsfeld.ord", col.names=cn)
freq <- data[,5]
data[,1] <- as.ordered(data[,1])
data[,2] <- as.ordered(data[,2])
data[,3] <- as.ordered(data[,3])
data[,4] <- as.ordered(data[,4])
vars <- data[,1:4]

nclass <- 2
nvar <- 4
nthresh <- 1

nameList <- names(vars)

class1 <- mxModel("Class1", 
            mxMatrix("Iden", name = "R", nrow = nvar, ncol = nvar, free=FALSE),
            mxMatrix("Full", name = "M", nrow = 1, ncol = nvar, free=FALSE),
            mxMatrix("Full", name = "ThresholdsClass1", nrow = 1, ncol = nvar, 
            	dimnames = list("Threshold",nameList), free=TRUE),
		  mxExpectationNormal(covariance="R", means="M", dimnames=nameList, 
				      thresholds="ThresholdsClass1"),
		  mxFitFunctionML(vector=TRUE))

class2 <- mxModel("Class2", 
            mxMatrix("Iden", name = "R", nrow = nvar, ncol = nvar, free=FALSE),
            mxMatrix("Full", name = "M", nrow = 1, ncol = nvar, free=FALSE),
            mxMatrix("Full", name = "ThresholdsClass2", nrow = 1, ncol = nvar, 
            	dimnames = list("Threshold",nameList), free=TRUE),
		  mxExpectationNormal(covariance="R", means="M", dimnames=nameList,
				      thresholds="ThresholdsClass2"),
		  mxFitFunctionML(vector=TRUE))

# Define the model

lcamodel <- mxModel("lcamodel", class1, class2, mxData(vars, type="raw"), 
            
# Create class membership probabilities, constrain them to be in the range 0-1, and make their sum equal 1.0
            mxMatrix("Full", name = "ClassMembershipProbabilities", nrow = nclass, ncol = 1, free=TRUE, 
			        labels = c(paste("pclass", 1:nclass, sep=""))),
            mxBounds(c(paste("pclass", 1:nclass, sep="")),0,1),
            mxConstraint(1 == sum(ClassMembershipProbabilities), name = 'sumConstraint'),

# Define the objective function
            mxAlgebra(-2*sum(freq * log(pclass1%x%Class1.objective + pclass2%x%Class2.objective)), 
            	name="lca"),
            mxFitFunctionAlgebra("lca")
)

# Run the job
model <- mxRun(lcamodel, suppressWarnings=TRUE)
summary(model)
model$matrices


# Check results against those hard-coded from old Mx:
negativeThresholds <- as.vector(mxEval(Class1.ThresholdsClass1, model))
positiveThresholds <- as.vector(mxEval(Class2.ThresholdsClass2, model))
if (negativeThresholds[[1]] > 0) {
	temp <- negativeThresholds
	negativeThresholds <- positiveThresholds
	positiveThresholds <- temp
}

omxCheckCloseEnough(negativeThresholds, c(-1.3247, -0.2909, -0.1466, -0.0044),.01)
omxCheckCloseEnough(positiveThresholds, c(0.1805, 0.9071, 1.3169, 1.5869),.01)
omxCheckCloseEnough(sort(as.vector(mxEval(ClassMembershipProbabilities, model))), c(0.4440, 0.5560),.01)
omxCheckCloseEnough(model$output$Minus2LogLikelihood, 4696.444, 0.01)
