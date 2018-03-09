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

library(OpenMx)

testmod1 <- mxModel(
	"NoJacobians",
	mxMatrix(type="Full",nrow=3,ncol=1,free=T,values=0.1,labels=paste("x",1:3,sep=""),lbound=0,name="X"),
	mxAlgebra( 3*X[1,1] + X[2,1] + X[3,1], name="fitfunc"),
	mxMatrix(type="Full",nrow=1,ncol=3,free=F,values=c(3,1,1),name="objgrad",dimnames=list(NULL,paste("x",1:3,sep=""))),
	mxConstraint(2*X[1,1] + X[2,1] + X[3,1] - 2 < 0,name="c1"),
	mxConstraint(X[1,1] - X[2,1] - X[3,1] + 1 < 0,name="c2"),
	mxFitFunctionAlgebra(algebra="fitfunc",gradient="objgrad")
)
testrun1 <- mxRun(testmod1)
summary(testrun1)
testrun1$fitfunction$result
testrun1$output$iterations
testrun1$output$evaluations


testmod2 <- mxModel(
	"Jacobians",
	mxMatrix(type="Full",nrow=3,ncol=1,free=T,values=0.1,labels=paste("x",1:3,sep=""),lbound=0,name="X"),
	mxAlgebra( 3*X[1,1] + X[2,1] + X[3,1], name="fitfunc"),
	mxMatrix(type="Full",nrow=1,ncol=3,free=F,values=c(3,1,1),name="objgrad",dimnames=list(NULL,paste("x",1:3,sep=""))),
	mxConstraint(2*X[1,1] + X[2,1] + X[3,1] - 2 < 0,name="c1",jac="jac1"),
	mxConstraint(X[1,1] - X[2,1] - X[3,1] + 1 < 0,name="c2",jac="jac2"),
	mxMatrix(type="Full",nrow=1,ncol=3,free=F,values=c(2,1,1),name="jac1",dimnames=list(NULL,paste("x",1:3,sep=""))),
	mxMatrix(type="Full",nrow=1,ncol=3,free=F,values=c(1,-1,-1),name="jac2",dimnames=list(NULL,paste("x",1:3,sep=""))),
	mxFitFunctionAlgebra(algebra="fitfunc",gradient="objgrad")
)
testrun2 <- mxRun(testmod2)
summary(testrun2)
testrun2$fitfunction$result
testrun2$output$iterations
testrun2$output$evaluations
omxCheckCloseEnough(mxEval(X,testrun1,T),mxEval(X,testrun2,T),1e-8)
omxCheckCloseEnough(mxEval(fitfunc,testrun1,T),mxEval(fitfunc,testrun2,T),1e-8)
omxCheckTrue(testrun1$output$evaluations > testrun2$output$evaluations)

