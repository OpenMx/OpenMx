#
#   Copyright 2007-2019 by the individuals mentioned in the source code history
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

# This script (written by Rob K.) is a regression test for a CSOLNP bug he repaired in April 2019 (2d34c67).
# It first run an unidentified exploratory factor analysis of unstandardized data, and then tries to 
# rotate the solution by minimizing an oblique rotation criterion, subject to the constraint
# that the model-expected covariance matrix remains the same as in the unrotated solution.
# The second MxModel (for rotation) is notable as an instance of the general case of a
# system of linearly dependent equality constraints lacking exact duplicates or proportional duplicates.
# CSOLNP is not expected to reach a sane solution; the test passes as long as CSOLNP
# does not segfault or raise C++ runtime errors.

library(OpenMx)
if(mxOption(NULL,"Default optimizer")!="CSOLNP"){stop("SKIP")}
set.seed(47402087)
data(HS.ability.data)

dataformx <- HS.ability.data[,7:30]
for(i in 1:24){
	dataformx[,i] <- as.double(dataformx[,i])
}
str(dataformx)
dataformx <- as.matrix(dataformx)
str(dataformx)

HSmodel <- mxModel(
	"Holzinger_and_Swineford_1939",
	mxData(observed=dataformx, type="raw"),
	mxMatrix(type="Full",nrow=1,ncol=24,free=T,values=colMeans(dataformx,na.rm=T),name="Mu"),
	mxMatrix(type="Full",nrow=24,ncol=3,free=T,values=0.1,name="Lambda"),
	mxMatrix(type="Stand",nrow=3,ncol=3,free=T,values=0.1,lbound=-0.9999,ubound=0.9999,name="Sxx"),
	mxMatrix(type="Diag",nrow=24,ncol=24,free=T,values=diag(cov(dataformx,use="pair")),lbound=0.0001,name="Psi2"),
	mxAlgebra(Lambda%&%Sxx + Psi2, name="xpecCov"),
	mxExpectationNormal(
		covariance="xpecCov",
		means="Mu",
		dimnames=colnames(dataformx)),
	mxFitFunctionML()
)
fitModel <- mxRun(HSmodel)
fitModel <- mxTryHard(fitModel,OKstatuscodes=c(0,1,5,6),exhaustive=T)

xpeccov1 <- mxEval(xpecCov,model=fitModel,compute=T)

rotmodel <- mxModel(
	"Holzinger_and_Swineford_1939",
	mxComputeSequence(steps=list(GD=mxComputeGradientDescent(engine=mxOption(NULL,"Default optimizer"),maxMajorIter=3000))),
	mxMatrix(type="Full",nrow=1,ncol=24,free=T,values=colMeans(dataformx,na.rm=T),name="Mu"),
	mxMatrix(type="Full",nrow=24,ncol=3,free=T,values=0.1,name="Lambda"),
	mxMatrix(type="Stand",nrow=3,ncol=3,free=T,values=0.1,lbound=-0.9999,ubound=0.9999,name="Sxx"),
	mxMatrix(type="Diag",nrow=24,ncol=24,free=T,values=diag(cov(dataformx,use="pair")),lbound=0.0001,name="Psi2"),
	mxMatrix(type="Symm",nrow=24,ncol=24,free=F,values=xpeccov1[lower.tri(xpeccov1,diag=T)],name="xpecCov1"),
	mxAlgebra(Lambda%&%Sxx + Psi2, name="xpecCov2"),
	mxConstraint(vech(xpecCov2) == vech(xpecCov1), name="myConstraint"),
	mxAlgebra(
		sum((Lambda[,1]%^%2)*(Lambda[,2]%^%2)) - 0.5/24*sum((Lambda[,1]%^%2))*sum((Lambda[,2]%^%2)),
		name="Part12"),
	mxAlgebra(
		sum((Lambda[,1]%^%2)*(Lambda[,3]%^%2)) - 0.5/24*sum((Lambda[,1]%^%2))*sum((Lambda[,3]%^%2)),
		name="Part13"),
	mxAlgebra(
		sum((Lambda[,2]%^%2)*(Lambda[,3]%^%2)) - 0.5/24*sum((Lambda[,2]%^%2))*sum((Lambda[,3]%^%2)),
		name="Part23"),
	mxAlgebra(Part12+Part13+Part23, name="fitfunc"),
	mxFitFunctionAlgebra(algebra="fitfunc",numObs=301,numStats=7224)
)
rotmodel <- omxSetParameters(model=rotmodel,labels=names(coef(fitModel)),values=coef(fitModel),free=c(rep(T,24),rep(T,99)))
rotfit <- mxRun(rotmodel)
rotfit <- mxModel(
	rotfit,
	mxComputeSequence(steps=list(GD=mxComputeGradientDescent(engine=mxOption(NULL,"Default optimizer"),maxMajorIter=3000)))
)
rotfit <- mxTryHard(rotfit,exhaustive=T)
