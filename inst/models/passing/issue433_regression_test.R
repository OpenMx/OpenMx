#
#   Copyright 2007-2026 by the individuals mentioned in the source code history
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
#   
#   
#   Regression test for issue #430, based on syntax contributed by Davood Tofighi.

library(OpenMx)
#SLSQP cannot reach the solution in one try, 
#and using `mxTryHard*()` would reach the solution on an attempt without any parameters necessarily starting at zero:
if(mxOption(NULL, 'Default optimizer') == "SLSQP") stop("SKIP")


#Generate data:
set.seed(1); n <- 100; x <- rnorm(n); m1 <- rnorm(n); y <- rnorm(n)
man <- c("x","m1","y")


mxOption(NULL,"Nudge zero starts","Yes")
mod <- mxModel(
	"c", type="RAM", manifestVars=man,
	mxPath("one", man, free=TRUE, values=0),
	mxPath("x","m1", free=TRUE, values=0, labels="b1"),
	mxPath("m1","y", free=TRUE, values=0, labels="b2"),
	mxPath(man, arrows=2, free=TRUE, values=1),
	mxData(data.frame(x,m1,y), type="raw"),
	mxAlgebra(b1*b2, name="ind"), mxConstraint(ind == 0, name="ind0", verbose=0))
omxCheckError(
	omxCheckWarning(mxRun(mod),""),
	"Expected `expression` to produce warnings."
)
fit <- mxRun(mod)
omxCheckTrue(abs(mxEval(ind,fit)[1,1]) < 1e-12)


mxOption(NULL,"Nudge zero starts","No")
mod2 <- mxModel(
	"c", type="RAM", manifestVars=man,
	mxPath("one", man, free=TRUE, values=0),
	mxPath("x","m1", free=TRUE, values=0, labels="b1"),
	mxPath("m1","y", free=TRUE, values=0, labels="b2"),
	mxPath(man, arrows=2, free=TRUE, values=1),
	mxData(data.frame(x,m1,y), type="raw"),
	mxAlgebra(b1*b2, name="ind"), mxConstraint(ind == 0, name="ind0", verbose=0))
fit2 <- omxCheckWarning(mxRun(mod2), "Ignoring constraint 'c.ind0[1]' because it does not depend on any free parameters")
omxCheckTrue(abs(mxEval(ind,fit2)[1,1]) > 1e-12)


mod3 <- mxModel(
	"c", type="RAM", manifestVars=man,
	mxPath("one", man, free=TRUE, values=0.1),
	mxPath("x","m1", free=TRUE, values=0.1, labels="b1"),
	mxPath("m1","y", free=TRUE, values=0.1, labels="b2"),
	mxPath(man, arrows=2, free=TRUE, values=1),
	mxData(data.frame(x,m1,y), type="raw"),
	mxAlgebra(b1*b2, name="ind"), mxConstraint(ind == 0, name="ind0", verbose=0))
omxCheckError(
	omxCheckWarning(mxRun(mod3),""),
	"Expected `expression` to produce warnings."
)
fit3 <- mxRun(mod3)
omxCheckTrue(abs(mxEval(ind,fit3)[1,1]) < 1e-12)




mxOption(reset=TRUE)
