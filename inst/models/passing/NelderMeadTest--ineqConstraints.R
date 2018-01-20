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
#Compare Nelder-Mead to the GD optimizer best at handling MxConstraints:
if(mxOption(NULL,"Default optimizer")!="SLSQP"){stop("SKIP")}
#The naive "soft" inequality method works reasonably well for inequalities:
foo <- mxComputeNelderMead(iniSimplexType="right", nudgeZeroStarts=FALSE, 
													 ineqConstraintMthd="soft", doPseudoHessian=T)
#foo$verbose <- 5L
plan <- omxDefaultComputePlan()
plan$steps <- list(foo,plan$steps$RE)

#Run with SLSQP:
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
testrun1$output$evaluations
testrun1$output$fit
testrun1$output$constraintFunctionValues

#Run with custom NM compute plan:
testmod2 <- mxModel(
	"NoJacobians",
	plan,
	mxMatrix(type="Full",nrow=3,ncol=1,free=T,values=c(0,0.5,0),labels=paste("x",1:3,sep=""),lbound=0,name="X"),
	mxAlgebra( 3*X[1,1] + X[2,1] + X[3,1], name="fitfunc"),
	mxMatrix(type="Full",nrow=1,ncol=3,free=F,values=c(3,1,1),name="objgrad",dimnames=list(NULL,paste("x",1:3,sep=""))),
	mxConstraint(2*X[1,1] + X[2,1] + X[3,1] - 2 < 0,name="c1"),
	mxConstraint(X[1,1] - X[2,1] - X[3,1] + 1 < 0,name="c2"),
	mxFitFunctionAlgebra(algebra="fitfunc",gradient="objgrad")
)
testrun2 <- mxRun(testmod2)
summary(testrun2)
testrun2$output$evaluations
testrun2$output$fit
#c2 is barely satisfied within feasibility tolerance:
mxEval(X[1,1] - X[2,1] - X[3,1] + 1, testrun2, T)

omxCheckCloseEnough(testrun2$output$fit+mxEval(X[1,1] - X[2,1] - X[3,1] + 1, testrun2, T), testrun1$output$fit, 1e-7)

#Backtracking:
foo <- mxComputeNelderMead(iniSimplexType="right", nudgeZeroStarts=FALSE, doPseudoHessian=T,
													 ineqConstraintMthd="eqMthd", eqConstraintMthd="backtrack")
#foo$verbose <- 5L
plan <- omxDefaultComputePlan()
plan$steps <- list(foo,plan$steps$RE)

testmod3 <- mxModel(
	"NoJacobians",
	plan,
	mxMatrix(type="Full",nrow=3,ncol=1,free=T,values=0.1,labels=paste("x",1:3,sep=""),lbound=0,name="X"),
	mxAlgebra( 3*X[1,1] + X[2,1] + X[3,1], name="fitfunc"),
	mxMatrix(type="Full",nrow=1,ncol=3,free=F,values=c(3,1,1),name="objgrad",dimnames=list(NULL,paste("x",1:3,sep=""))),
	mxConstraint(2*X[1,1] + X[2,1] + X[3,1] - 2 < 0,name="c1"),
	mxConstraint(X[1,1] - X[2,1] - X[3,1] + 1 < 0,name="c2"),
	mxFitFunctionAlgebra(algebra="fitfunc",gradient="objgrad")
)
mxEval(X[1,1] - X[2,1] - X[3,1] + 1, testmod3, T) #<--More than feas tolerance
testrun3 <- mxRun(testmod3)
summary(testrun3)
testrun3$output$evaluations
testrun3$output$fit
#c2 is barely satisfied within feasibility tolerance:
mxEval(X[1,1] - X[2,1] - X[3,1] + 1, testrun3, T)

omxCheckCloseEnough(testrun3$output$fit+mxEval(X[1,1] - X[2,1] - X[3,1] + 1, testrun3, T), testrun1$output$fit, 1e-7)


#l1p:
foo <- mxComputeNelderMead(iniSimplexType="right", nudgeZeroStarts=FALSE, doPseudoHessian=T,
													 ineqConstraintMthd="eqMthd", eqConstraintMthd="l1p")
#foo$verbose <- 5L
plan <- omxDefaultComputePlan()
plan$steps <- list(foo,plan$steps$RE)

testmod4 <- mxModel(
	"NoJacobians",
	plan,
	mxMatrix(type="Full",nrow=3,ncol=1,free=T,values=0.1,labels=paste("x",1:3,sep=""),lbound=0,name="X"),
	mxAlgebra( 3*X[1,1] + X[2,1] + X[3,1], name="fitfunc"),
	mxMatrix(type="Full",nrow=1,ncol=3,free=F,values=c(3,1,1),name="objgrad",dimnames=list(NULL,paste("x",1:3,sep=""))),
	mxConstraint(2*X[1,1] + X[2,1] + X[3,1] - 2 < 0,name="c1"),
	mxConstraint(X[1,1] - X[2,1] - X[3,1] + 1 < 0,name="c2"),
	mxFitFunctionAlgebra(algebra="fitfunc",gradient="objgrad")
)
testrun4 <- mxRun(testmod4)
summary(testrun4)
omxCheckCloseEnough(testrun4$output$fit+mxEval(X[1,1] - X[2,1] - X[3,1] + 1, testrun4, T), testrun1$output$fit, 1e-7)
#PseudoHessian shouldn't be calculated when using l1p:
omxCheckTrue(is.null(testrun4$compute$steps[[1]]$output$pseudoHessian))


#GDsearch:
foo <- mxComputeNelderMead(iniSimplexType="right", nudgeZeroStarts=FALSE, doPseudoHessian=T,
													 ineqConstraintMthd="eqMthd", eqConstraintMthd="GDsearch")
#foo$verbose <- 5L
plan <- omxDefaultComputePlan()
plan$steps <- list(foo,plan$steps$RE)

testmod5 <- mxModel(
	"NoJacobians",
	plan,
	mxMatrix(type="Full",nrow=3,ncol=1,free=T,values=0.1,labels=paste("x",1:3,sep=""),lbound=0,name="X"),
	mxAlgebra( 3*X[1,1] + X[2,1] + X[3,1], name="fitfunc"),
	mxMatrix(type="Full",nrow=1,ncol=3,free=F,values=c(3,1,1),name="objgrad",dimnames=list(NULL,paste("x",1:3,sep=""))),
	mxConstraint(2*X[1,1] + X[2,1] + X[3,1] - 2 < 0,name="c1"),
	mxConstraint(X[1,1] - X[2,1] - X[3,1] + 1 < 0,name="c2"),
	mxFitFunctionAlgebra(algebra="fitfunc",gradient="objgrad")
)
testrun5 <- mxRun(testmod5)
summary(testrun5)
omxCheckCloseEnough(testrun5$output$fit+mxEval(X[1,1] - X[2,1] - X[3,1] + 1, testrun5, T), testrun1$output$fit, 1e-7)
