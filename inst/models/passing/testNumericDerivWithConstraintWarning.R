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

set.seed(1611150)
x <- matrix(rnorm(1000,sd=2))
colnames(x) <- "x"

plan <- omxDefaultComputePlan(modelName="mod")

varmod <- mxModel(
	"mod",
	plan,
	mxData(observed=x,type="raw"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0,labels="mu",name="Mu"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=4,labels="sigma2",name="Sigma2",lbound=0),
	mxExpectationNormal(covariance="Sigma2",means="Mu",dimnames=c("x")),
	mxAlgebra(sqrt(Sigma2),name="Sigma"),
	mxFitFunctionML(),
	mxConstraint(Sigma2 - Mu > 0, name="pointless")
)
#NPSOL also warns about status code 1:
if(mxOption(NULL,"Default optimizer")!='NPSOL'){
	omxCheckWarning(varrun <- mxRun(varmod), 
									"due to presence of MxConstraints, Hessian matrix and standard errors may not be valid for statistical-inferential purposes")
#(^^^Ironically, the SEs are probably fine in this case, since the only constraint is an inactive inequality...)
	summary(varrun)
}


plan2 <- omxDefaultComputePlan(modelName="mod")
plan2$steps$GD <- mxComputeNelderMead(iniSimplexType="right", nudgeZeroStarts=FALSE, 
																			ineqConstraintMthd="soft", doPseudoHessian=T)
varmod2 <- mxModel(
	"mod",
	plan2,
	mxData(observed=x,type="raw"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0,labels="mu",name="Mu"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=4,labels="sigma2",name="Sigma2",lbound=0),
	mxExpectationNormal(covariance="Sigma2",means="Mu",dimnames=c("x")),
	mxAlgebra(sqrt(Sigma2),name="Sigma"),
	mxFitFunctionML(),
	mxConstraint(Sigma2 - Mu > 0, name="pointless")
)
varrun2 <- mxRun(varmod2)
varrun2$output$standardErrors
sqrt(diag(chol2inv(chol(varrun2$compute$steps$GD$output$pseudoHessian))))
