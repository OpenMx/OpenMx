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

# Created by Mike Cheung

require(OpenMx)

yi <- c(-0.264,-0.230,0.166,0.173,0.225,0.291,0.309,0.435,0.476,0.617,0.651,0.718,0.740,0.745,0.758,0.922,0.938,0.962,1.522,1.844)
vi <- c(0.086,0.106,0.055,0.084,0.071,0.078,0.051,0.093,0.149,0.095,0.110,0.054,0.081,0.084,0.087,0.103,0.113,0.083,0.100,0.141)
my.df <- cbind(yi,vi)
test <- mxModel("test", type="default",
	mxMatrix("Zero", ncol=1, nrow=1, free=F, name="Amat"),
		mxAlgebra(Amat[1,1], name="A"),  # just to test with A as an algebra
	mxMatrix("Full", ncol=1, nrow=1, free=F, values=0,   labels="data.vi", name="V"),
	mxMatrix("Full", ncol=1, nrow=1, free=T, values=0.1, lbound=0.0000001, name="Tau"),
	mxMatrix("Full", ncol=1, nrow=1, free=T, values=0,   name="M"),
	mxMatrix("Iden", ncol=1, nrow=1, name="F"),
	mxAlgebra(V+Tau, name="S"),
	mxFitFunctionML(),mxExpectationRAM("A", "S", "F", "M", dimnames=c("yi")),
	mxData(observed=my.df, type="raw")
)
out <- mxRun(test, suppressWarnings=TRUE)
omxCheckCloseEnough(mxEval(objective, out), 27.8, 0.01)
