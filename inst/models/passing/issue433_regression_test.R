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
#   Regression test for issue #433, based on syntax contributed by Davood Tofighi.

library(OpenMx)
#SLSQP cannot reach the solution in one try, 
#and using `mxTryHard*()` would reach the solution on an attempt without any parameters necessarily starting at zero;
#NPSOL emits extra warnings about status GREEN / status code 1:
if(mxOption(NULL, 'Default optimizer') != "CSOLNP") stop("SKIP")


#Generate data:
set.seed(1); n <- 100; x <- rnorm(n); m1 <- rnorm(n); y <- rnorm(n)
man <- c("x","m1","y")


#Test non-default value of this mxOption:
mxOption(NULL,"Drop extraneous MxConstraints","No")
mod <- mxModel(
	"c", type="RAM", manifestVars=man,
	mxPath("one", man, free=TRUE, values=0),
	mxPath("x","m1", free=TRUE, values=0, labels="b1"),
	mxPath("m1","y", free=TRUE, values=0, labels="b2"),
	mxPath(man, arrows=2, free=TRUE, values=1),
	mxData(data.frame(x,m1,y), type="raw"),
	mxAlgebra(b1*b2, name="ind"), mxConstraint(ind == 0, name="ind0", verbose=0))
#"ind0" is not spuriously dropped:
omxCheckError(
	omxCheckWarning(mxRun(mod),""),
	"Expected `expression` to produce warnings."
)
fit <- mxRun(mod)
#Constraint is respected:
omxCheckTrue(abs(mxEval(ind,fit)[1,1]) < 1e-12)


#Restore on-load default value of this mxOption:
mxOption(NULL,"Drop extraneous MxConstraints","Yes")
mod2 <- mxModel(
	"c", type="RAM", manifestVars=man,
	mxPath("one", man, free=TRUE, values=0),
	mxPath("x","m1", free=TRUE, values=0, labels="b1"),
	mxPath("m1","y", free=TRUE, values=0, labels="b2"),
	mxPath(man, arrows=2, free=TRUE, values=1),
	mxData(data.frame(x,m1,y), type="raw"),
	mxAlgebra(b1*b2, name="ind"), mxConstraint(ind == 0, name="ind0", verbose=0))
#OpenMx drops the spuriously extraneous constraint, and warns:
fit2 <- omxCheckWarning(mxRun(mod2), "Ignoring constraint 'c.ind0[1]' because it does not depend on any free parameters")
#Constraint is not respected:
omxCheckTrue(abs(mxEval(ind,fit2)[1,1]) > 1e-12)


#Use nonzero start values:
mod3 <- mxModel(
	"c", type="RAM", manifestVars=man,
	mxPath("one", man, free=TRUE, values=0.1),
	mxPath("x","m1", free=TRUE, values=0.1, labels="b1"),
	mxPath("m1","y", free=TRUE, values=0.1, labels="b2"),
	mxPath(man, arrows=2, free=TRUE, values=1),
	mxData(data.frame(x,m1,y), type="raw"),
	mxAlgebra(b1*b2, name="ind"), mxConstraint(ind == 0, name="ind0", verbose=0))
#"ind0" is not dropped:
omxCheckError(
	omxCheckWarning(mxRun(mod3),""),
	"Expected `expression` to produce warnings."
)
fit3 <- mxRun(mod3)
#Constraint is respected:
omxCheckTrue(abs(mxEval(ind,fit3)[1,1]) < 1e-12)


#Jiggling the start values keeps "ind0" from being dropped:
omxCheckError(
	omxCheckWarning(mxRun(mxJiggle(model=mod2,classic=T)),""),
	"Expected `expression` to produce warnings."
)
fit2a <- mxRun(mxJiggle(model=mod2,classic=T))
omxCheckTrue(abs(mxEval(ind,fit2a)[1,1]) < 1e-12)
set.seed(1); omxCheckError(
	omxCheckWarning(mxRun(mxJiggle(model=mod2,classic=F)),""),
	"Expected `expression` to produce warnings."
)
set.seed(1); fit2b <- mxRun(mxJiggle(model=mod2,classic=F))
omxCheckTrue(abs(mxEval(ind,fit2b)[1,1]) < 1e-12)


if(0){
	mxOption(key="Feasibility tolerance",value=1e-5)
	set.seed(1); fit2c <- mxTryHard(mod2)
	omxCheckTrue(abs(mxEval(ind,fit2c)[1,1]) < 1e-12)
}




mxOption(reset=TRUE)
