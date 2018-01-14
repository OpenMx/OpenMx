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
## Create n=100 binary data with p=.6
my.df <- data.frame(x=c(rep(0,each=40),rep(1,each=60)))
## Starting value for the threshold
st <- qnorm(colMeans(my.df))
my.df$x <- mxFactor(my.df$x, levels=c(0,1))

prop.Model <- mxModel("Binary variable", mxData(my.df, type="raw"),
mxMatrix(type="Full", nrow=1, ncol=2, values=c(0,1), free=FALSE, name="A"),
mxMatrix(type="Diag", nrow=2, ncol=2, values=c(1,1), free=FALSE, name="S"),
mxMatrix(type="Full", nrow=1, ncol=2, values=c(1,0), free=FALSE, name="F"),
mxMatrix(type="Full", nrow=1, ncol=2, values=c(0,0), free=FALSE, name="M"),
mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, values=st, labels="th1", name="thresh"),
mxFitFunctionML(),mxExpectationRAM("A", "S", "F", "M", dimnames=c("x","f1"),
thresholds="thresh", threshnames=c("x"))
)
omxCheckError(mxRun(prop.Model),
              "RAM matrices 'Binary variable.S' and 'Binary variable.A' must have the same dimensions")

noData <- mxModel("No data",
		  mxMatrix(type="Full", nrow=2, ncol=2, values=c(0,1), free=FALSE, name="A"),
		  mxMatrix(type="Diag", nrow=2, ncol=2, values=c(1,1), free=FALSE, name="S"),
		  mxMatrix(type="Full", nrow=1, ncol=2, values=c(1,0), free=FALSE, name="F"),
		  mxMatrix(type="Full", nrow=1, ncol=2, values=c(0,0), free=FALSE, name="M"),
		  mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, values=st, labels="th1", name="thresh"),
		  mxFitFunctionML(),
		  mxExpectationRAM(M="M", dimnames=c("x","f1"),
				   thresholds="thresh", threshnames=c("x"))
)
omxCheckError(mxRun(noData), paste("The RAM expectation function does not",
				   "have a dataset associated with it in model 'No data'"))
