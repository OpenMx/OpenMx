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


#This script is used to test the definition variable functionality in OpenMx
#The definition variable in this example is dichotomous, and describes two different groups
#These two groups are measured on two variables, x and y
#The group with a definition value of 1 has a perfect relationship between x and y
#The group with a definition value of 0 has no relationship between x and y (beyond chance)
#The definition variable will then be used to define the covarinace between x and y


#make some data!
x1<-rnorm(50)
y1<-x1+rnorm(50, sd=1)

x2<-rnorm(50)
y2<-rnorm(50)+rnorm(50, sd=1)

#put them both together, add a definition variable, and make an MxData object
x<-c(x1,x2)
y<-c(y1,y2)
def<-rep(c(.88,0),each=50)
data<-mxData(as.matrix(data.frame(x,y,def)), type="raw")

#define the model: we'll just use an S matrix and let A and F drop out
#as currently specified, this would fit a zero df model to a 2x2 covariance matrix
S <- mxMatrix("Symm", values=c(2, 0.5, 3), labels=c(NA, "data.def", NA),
	lbound=c(0.001, NA, 0.001), free=TRUE, nrow=2, ncol=2, name="S",
	dimnames = list(c('x','y'), c('x','y')))

M <- mxMatrix("Zero", nrow = 1, ncol = 2, name = "M", dimnames = list(NULL, c('x','y'))) 

#define the model, including a FIML objective function, which will optimize the matrix S
model<-mxModel("model", mxFIMLObjective("S", "M"), data, S, M)

#run the model
run<-mxRun(model)

observed <- diag(cov(cbind(x,y)))
computed <- c(run@output$estimate)

omxCheckCloseEnough(observed, computed, 0.1)


