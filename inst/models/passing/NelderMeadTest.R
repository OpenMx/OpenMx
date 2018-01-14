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
#This test involves CIs, so use it with SLSQP, since the inequality-constrained formulation is theoretically more correct:
if(mxOption(NULL,"Default optimizer")!="SLSQP"){stop("SKIP")}
#Need to use stricter convergence tolerances to avoid status Red:
foo <- mxComputeNelderMead(iniSimplexType="smartRight", nudgeZeroStarts=FALSE, xTolProx=1e-8, fTolProx=1e-8,
													 doPseudoHessian=T)
#foo$verbose <- 5L
plan <- omxDefaultComputePlan(intervals=T)
plan$steps$GD <- foo
plan$steps$CI$plan <- mxComputeNelderMead()
plan$steps$CI$constraintType <- "none"

#Simulate data:
set.seed(1611150)
x <- matrix(rnorm(1000,sd=2))
colnames(x) <- "x"

#Summary statistics:
print(mean(x))
print(var(x))

#Run with SLSQP:
varmodGD <- mxModel(
	"mod",
	mxData(observed=x,type="raw"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0,labels="mu",name="Mu"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=4,labels="sigma2",name="Sigma2",lbound=0),
	mxExpectationNormal(covariance="Sigma2",means="Mu",dimnames=c("x")),
	mxAlgebra(sqrt(Sigma2),name="Sigma"),
	mxCI(c("mu","sigma2")),
	mxFitFunctionML()
)
varrunGD <- mxRun(varmodGD,intervals=T)

#Run with custom NM compute plan:
varmod <- mxModel(
	"mod",
	plan,
	mxData(observed=x,type="raw"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0,labels="mu",name="Mu"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=4,labels="sigma2",name="Sigma2",lbound=0),
	mxExpectationNormal(covariance="Sigma2",means="Mu",dimnames=c("x")),
	mxAlgebra(sqrt(Sigma2),name="Sigma"),
	mxCI(c("mu","sigma2")),
	mxFitFunctionML()
)
varrun <- mxRun(varmod,intervals=T)

#Tests:
c(mean(x),var(x)*999/1000)
varrunGD$output$estimate
varrun$output$estimate
omxCheckCloseEnough(varrun$output$estimate, c(mean(x),var(x)*999/1000), 1e-5)

varrunGD$output$fit
varrun$output$fit
omxCheckCloseEnough(varrunGD$output$fit, varrun$output$fit, 1e-4)

varrunGD$output$standardErrors
varrun$output$standardErrors
omxCheckCloseEnough(varrunGD$output$standardErrors, varrun$output$standardErrors, 1e-5)
omxCheckCloseEnough(
	sqrt(diag(chol2inv(chol(varrun$compute$steps[[1]]$output$pseudoHessian)))),
	as.vector(varrunGD$output$standardErrors),
	1e-04)

varrunGD$output$confidenceIntervals
varrun$output$confidenceIntervals
omxCheckCloseEnough(
	varrunGD$output$confidenceIntervals,
	varrun$output$confidenceIntervals,
	0.01
)

omxCheckTrue(length(varrun$compute$steps$GD$output$paramNames))
omxCheckTrue(length(varrun$compute$steps$GD$output$finalSimplexMat))
omxCheckTrue(length(varrun$compute$steps$GD$output$finalFitValues))
omxCheckTrue(length(varrun$compute$steps$GD$output$finalVertexInfeas))
omxCheckTrue(length(varrun$compute$steps$GD$output$pseudoHessian))
omxCheckTrue(length(varrun$compute$steps$GD$output$simplexGradient))
omxCheckTrue(length(varrun$compute$steps$GD$output$rangeProximityMeasure))
omxCheckTrue(length(varrun$compute$steps$GD$output$domainProximityMeasure))
omxCheckTrue(length(varrun$compute$steps$GD$output$penalizedFit))

#Try using inequality-constrained formulation of CI problem:
plan$steps$CI$constraintType <- "ineq"
varmod2 <- mxModel(
	"mod",
	plan,
	mxData(observed=x,type="raw"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0,labels="mu",name="Mu"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=4,labels="sigma2",name="Sigma2",lbound=0),
	mxExpectationNormal(covariance="Sigma2",means="Mu",dimnames=c("x")),
	mxAlgebra(sqrt(Sigma2),name="Sigma"),
	mxCI(c("mu","sigma2")),
	mxFitFunctionML()
)
varrun2 <- mxRun(varmod2,intervals=T)
varrunGD$output$confidenceIntervals
varrun2$output$confidenceIntervals
omxCheckCloseEnough(
	varrunGD$output$confidenceIntervals,
	varrun2$output$confidenceIntervals,
	0.01 
)
#The change in fit is off by 0.05, which is the on-load default feasibility tolerance:
omxCheckCloseEnough(varrun2$compute$steps$CI$output$detail$fit - varrun2$output$fit - 0.05, rep(qchisq(0.95,1),4), 1e-6)
