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

varmod <- mxModel(
	"mod",
	mxData(observed=x,type="raw"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0,labels="mu",name="Mu"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=4,labels="sigma2",name="Sigma2",lbound=0),
	mxExpectationNormal(covariance="Sigma2",means="Mu",dimnames=c("x")),
	mxAlgebra(sqrt(Sigma2),name="Sigma"),
	mxFitFunctionML()
)
varrun <- mxRun(varmod)

sdmod <- mxModel(
	"mod",
	mxData(observed=x,type="raw"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0,labels="mu",name="Mu"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=2,labels="sigma",name="Sigma",lbound=0),
	mxExpectationNormal(covariance="Sigma2",means="Mu",dimnames=c("x")),
	mxAlgebra(Sigma^2,name="Sigma2"),
	mxFitFunctionML()
)
sdrun <- mxRun(sdmod)

omxCheckCloseEnough(
	varrun$output$standardErrors[2],
	mxSE(x=Sigma2,model=sdrun),
	1e-6
)

omxCheckCloseEnough(
	sdrun$output$standardErrors[2],
	mxSE(x=Sigma,model=varrun),
	1e-6
)

omxCheckCloseEnough(mxSE(sigma^2, sdrun), mxSE(sigma2, varrun), 1e-6)

