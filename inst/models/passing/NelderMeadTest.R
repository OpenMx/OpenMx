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
foo <- mxComputeNelderMead()
foo$verbose <- 5L
plan <- omxDefaultComputePlan()
plan$steps <- list(foo,plan$steps$RE)

set.seed(1611150)
x <- matrix(rnorm(1000,sd=2))
colnames(x) <- "x"

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
omxCheckError(mxRun(varmod),"NelderMeadOptimizerContext::invokeNelderMead() : so far, so good")
