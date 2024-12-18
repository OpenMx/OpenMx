#
#   Copyright 2007-2024 by the individuals mentioned in the source code history
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

# Adapted from script by Timothy C. Bates

require(OpenMx)
# This script doesn't test anything optimization-related, so it doesn't need 
# to be run with anything other than the on-load default:
if(mxOption(NULL, 'Default optimizer') != "SLSQP") stop("SKIP")
data(myGrowthMixtureData)

class1 <- mxModel(
	"Class1", type="RAM", manifestVars=c("x1","x2","x3","x4","x5"), latentVars=c("intercept","slope"),
	mxPath(from=c("x1","x2","x3","x4","x5"), arrows=2, free=TRUE, values = c(1,1,1,1,1),labels=c("residual","residual","residual","residual","residual") ),
	mxPath(from=c("intercept","slope"), arrows=2, connect="unique.pairs", free=TRUE, values=c(1,.4,1), labels=c("vari1","cov1","vars1") ),
	mxPath(from="intercept", to=c("x1","x2","x3","x4","x5"), arrows=1, free=FALSE, values=c(1,1,1,1,1) ),
	mxPath(from="slope", to=c("x1","x2","x3","x4","x5"), arrows=1, free=FALSE, values=c(0,1,2,3,4) ),
	mxPath(from="one", to=c("x1","x2", "x3", "x4","x5"), arrows=1, free=FALSE, values=c(0,0,0,0,0) ),
	mxPath(from="one", to=c("intercept","slope"), arrows=1, free=TRUE,  values=c(0,-1), labels=c("meani1","means1") ),
	mxFitFunctionML(vector=TRUE)
)

class2 <- mxModel(
	class1, name="Class2",
	mxPath(from=c("intercept","slope"), arrows=2, connect="unique.pairs", free=TRUE, values=c(1,.5,1), labels=c("vari2","cov2","vars2") ),
	mxPath(from="one", to=c("intercept", "slope"), arrows=1, free=TRUE, values=c(5,1), labels=c("meani2","means2") )
)

m1 <- mxModel(
	"Growth Mixture Model", mxData(myGrowthMixtureData, type="raw" ), class1, class2,
	mxMatrix(name="Props", type="Full", nrow=2, ncol=1, free=c(TRUE, FALSE), values=1, lbound=0.001, labels = c("p1","p2")),
	mxExpectationMixture(components=c('Class1', 'Class2'), weights='Props', scale='sum'),
	mxFitFunctionML()
)

old <- names(omxGetParameters(m1, indep = TRUE, free = TRUE))
# Should work without error:
m2  <- omxSetParameters(m1, old, free = FALSE, strict = TRUE, name = "_fixed")
# Should run without error:
m2 <- mxRun(m2)
