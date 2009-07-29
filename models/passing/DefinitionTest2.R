#
#   Copyright 2007-2009 The OpenMx Project
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

#Definition Variable Test 1
#Author: Ryne Estabrook
#Date: 12 May 2009
#Modified: Mike Neale
#Date: July 28 2009


#This script is used to test the definition variable functionality in OpenMx
#The definition variable in this example is dichotomous, and describes two different groups
#These two groups are measured on two variables, x and y
#The group with a definition value of 1 has a covariance of 1 between x and y
#The group with a definition value of 0 has no relationship between x and y (beyond chance)
#The definition variable is then used to define the covariance between x and y


#make some data!
x1<-rnorm(500)
y1<-x1+rnorm(500, sd=1)

x2<-rnorm(500)
y2<-rnorm(500)+rnorm(500, sd=1)

#put them both together, add a definition variable, and make an MxData object
x<-c(x1,x2)
y<-c(y1,y2)
def<-rep(c(1,0),each=500)
selvars<-c("x","y")


# Three covariance model matrices: 
#  "cov" for the zero relationship group
#  "def" for the definition variable, 
#   and "beta" for estimating difference between groups' covariances
# One common mean vector, "M"

model<-mxModel("model", mxFIMLObjective("S", "M"), 
				mxData(as.matrix(data.frame(x,y,def)), type="raw"),
				mxMatrix("Symm", nrow=2, ncol=2, free=FALSE, values=c(0, 0, 0), labels=c(NA, "data.def", NA),
					dimnames=list(selvars,selvars), name="def"),
				mxMatrix("Symm", nrow=2, ncol=2, free=c(FALSE,TRUE,TRUE,FALSE), values=c(0, 0.0001, .0001, 0),
#					lbound=c(NA, -1, -1, NA), ubound=c(NA, 1, 1, NA),
					dimnames=list(selvars,selvars), name="beta"),
				mxMatrix("Symm", nrow=2, ncol=2, free=TRUE, values=c(1, 0, 1),
					lbound=c(0.001, NA, 0.001), dimnames=list(selvars,selvars), name="cov"),
				mxMatrix("Full", nrow = 1, ncol = 2, free=TRUE, dimnames=list(NULL,selvars), name = "M"),
				mxAlgebra(cov+beta*def, name="S", dimnames=list(selvars,selvars))
			)
#define the model, including a FIML objective function, which will optimize the matrix S

#run the model
run<-mxRun(model)

#check results
ObsCov1<-(cov(cbind(x1,y1)))
ObsCov2<-(cov(cbind(x2,y2)))

#Note: this check may fail because due to laziness we use variance of group 2 only 
#      whereas variances estimated by the model are for the common group only
observed <- c(cov(x1,y1),ObsCov2[lower.tri(ObsCov2,diag=TRUE)],mean(x),mean(y))
estimated <- c(run@output$estimate)

omxCheckCloseEnough(observed, estimated, 0.1)


