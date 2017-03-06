#
#   Copyright 2007-2017 The OpenMx Project
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
#No need to run this test with other than the on-load default GD optimizer:
if(mxOption(NULL,"Default optimizer")!="CSOLNP"){stop("SKIP")}
#Need to use stricter convergence tolerances to avoid status Red:
foo <- mxComputeNelderMead(iniSimplexType="smartRight", nudgeZeroStarts=FALSE, xTolProx=1e-8, fTolProx=1e-8,
													 doPseudoHessian=T)
#foo$verbose <- 5L
plan <- omxDefaultComputePlan()
plan$steps$GD <- foo

#Simulate data:
set.seed(1611150)
x <- matrix(rnorm(1000,sd=2))
colnames(x) <- "x"

#Summary statistics:
print(mean(x))
print(var(x))

#Run with CSOLNP:
varmodGD <- mxModel(
	"mod",
	mxData(observed=x,type="raw"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0,labels="mu",name="Mu"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=4,labels="sigma2",name="Sigma2",lbound=0),
	mxExpectationNormal(covariance="Sigma2",means="Mu",dimnames=c("x")),
	mxAlgebra(sqrt(Sigma2),name="Sigma"),
	mxFitFunctionML()
)
varrunGD <- mxRun(varmodGD)

#Run with custom NM compute plan:
varmod <- mxModel(
	"mod",
	plan,
	mxData(observed=x,type="raw"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0,labels="mu",name="Mu"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=4,labels="sigma2",name="Sigma2",lbound=0),
	mxExpectationNormal(covariance="Sigma2",means="Mu",dimnames=c("x")),
	mxAlgebra(sqrt(Sigma2),name="Sigma"),
	mxFitFunctionML()
)
varrun <- mxRun(varmod)

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
