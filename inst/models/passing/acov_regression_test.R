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

library(OpenMx)

S <- matrix(
	c(2,1,1,1,2,1,1,1,2),
	ncol=3,nrow=3,
	dimnames=list(c("V1","V2","V3"),c("V1","V2","V3"))
)
W <- diag(6)

mxdataobj <- mxData(S, numObs = 100, means = NA,
                    type = "acov",acov=W, fullWeight=W)

factormod <- mxModel(
	"1factor",
	mxdataobj,
	mxMatrix(type="Full",nrow=3,ncol=1,free=T,values=0.1,name="Lambda"),
	mxMatrix(type="Stand",nrow=1,ncol=1,free=F,values=1,lbound=-1,ubound=1,name="Rxx"),
	mxMatrix(type="Diag",nrow=3,free=T,values=diag(S),lbound=0.0001,name="Psi2"),
	mxAlgebra( (Lambda%*%Rxx%*%t(Lambda)) + Psi2, name="xpecCov"),
	mxAlgebra(sqrt(diag2vec(xpecCov)), name="SDs"),
	mxAlgebra( 
		Lambda / SDs, name="standardizedLoadings"),
	mxExpectationNormal(
		covariance="xpecCov",
		dimnames=colnames(S)),
	mxFitFunctionWLS(type="WLS")
)
factorfit <- mxRun(factormod)

summary(factorfit)

omxCheckCloseEnough(as.vector(mxEval(Lambda,factorfit,T)),c(1,1,1),1e-6)
omxCheckCloseEnough(as.vector(mxEval(diag2vec(Psi2),factorfit,T)),c(1,1,1),1e-3)

#---

omxCheckError(mxRun(mxModel(factormod, mxData(S, numObs = 100, means = NA,
	type = "acov",acov=W))),
	"MxComputeStandardError: terribly sorry, master, but '1factor.data' does not include the full weight matrix hence standard errors cannot be computed")
