#
#   Copyright 2007-2020 by the individuals mentioned in the source code history
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

#The quadratic fitfunction in this model has a solution on a parameter boundary, a zero-gradient point outside the feasible space, and a
#Hessian matrix that is nowhere PD.

library(OpenMx)
library(testthat)

mxOption(key="feasibility tolerance", value = .00001)

startvals <- c(-5.1, 2.9)

plan <- omxDefaultComputePlan()
plan$steps <- list(plan$steps$GD)
#plan$steps[[1]]$verbose <- 1L

m1 <- mxModel(
	"BukinN2",
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=startvals[1],labels="x1",lbound=-15,ubound=-5,name="X1"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=startvals[2],labels="x2",lbound=-3,ubound=3,name="X2"),
	mxAlgebra( 100*( (X2^2) - 0.01*(X1^2) + 1) + 0.01*(X1+10)^2,
						 name="BukinN2Func"),
	#Interestingly, given an analytic gradient to NPSOL and SLSQP makes them FAIL to find the correct solution...
	mxAlgebra(cbind(-1.98*X1 + 0.2, 200*X2), name="grad", dimnames=list(NULL,c("x1","x2"))),
	mxFitFunctionAlgebra(algebra="BukinN2Func")#,gradient="grad")
)
m1 <- mxRun(m1)
summary(m1)
m1$output$gradient
omxCheckCloseEnough(coef(m1), c(-15,0), 0.1)
omxCheckCloseEnough(m1$output$fit, -124.7500, 5e-5)
omxCheckCloseEnough(m1$output$gradient[1],29.9,0.01)
omxCheckTrue(m1$output$status$code %in% c(0,1))


m2 <- mxModel(
	"BukinN2",
	plan,
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=startvals[1],labels="x1",lbound=-15,ubound=-5,name="X1"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=startvals[2],labels="x2",lbound=-3,ubound=3,name="X2"),
	mxAlgebra( 100*( (X2^2) - 0.01*(X1^2) + 1) + 0.01*(X1+10)^2,
						 name="BukinN2Func"),
	#Interestingly, given an analytic gradient to NPSOL and SLSQP makes them FAIL to find the correct solution...
	mxAlgebra(cbind(-1.98*X1 + 0.2, 200*X2), name="grad", dimnames=list(NULL,c("x1","x2"))),
	mxFitFunctionAlgebra(algebra="BukinN2Func")#,gradient="grad")
)
m2 <- mxRun(m2)
summary(m2)
m2$output$gradient
omxCheckCloseEnough(coef(m2), c(-15,0), 0.1)
omxCheckCloseEnough(m2$output$fit, -124.7500, 5e-5)
omxCheckTrue(m2$output$status$code %in% c(0,1))


#All 3 optimizers appear robust to redundant INequality constraints, even with infeasible start values:
startvals <- c(-4.9, 3.1)
m3 <- mxModel(
	"BukinN2",
	plan,
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=startvals[1],labels="x1",name="X1"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=startvals[2],labels="x2",name="X2"),
	mxConstraint(X1 > -15, name="l1"),
	mxConstraint(2*X1 > -30, name="redundant"),
	mxConstraint(X1 < -5, name="u1"),
	mxConstraint(X2 > -3, name="l2"),
	mxConstraint(X2 < 3, name="u2"),
	mxAlgebra( 100*( (X2^2) - 0.01*(X1^2) + 1) + 0.01*(X1+10)^2,
						 name="BukinN2Func"),
	mxAlgebra(cbind(-1.98*X1 + 0.2, 200*X2), name="grad", dimnames=list(NULL,c("x1","X2"))),
	mxFitFunctionAlgebra(algebra="BukinN2Func",gradient="grad")
)
expect_error(mxRun(m3), "'BukinN2.fitfunction' has a derivative entry for unrecognized parameter 'X2'")
colnames(m3$grad) <- c("x1","x2")
m3 <- mxRun(m3)
m3 <- mxRun(m3)
summary(m3)
omxCheckCloseEnough(coef(m3), c(-15,0), 0.1)
omxCheckCloseEnough(m3$output$fit, -124.7501, .0002)

m4 <- omxSetParameters(model=m1,labels="x2",free=F,values=0)
m4 <- mxRun(m4)
# Should still be TRUE if the sole free parameter has an active bound:
omxCheckTrue(m4$output$infoDefinite)

m5 <- mxModel(
	"BukinN2",
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=startvals[1],labels="x1",lbound=-15,ubound=-5,name="X1"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=startvals[2],labels="x2",lbound=-3,ubound=3,name="X2"),
	mxAlgebra( 100*( (X2^2) - 0.01*(X1^2) + 1) + 0.01*(X1+10)^2,
						 name="BukinN2Func"),
	#Interestingly, given an analytic gradient to NPSOL and SLSQP makes them FAIL to find the correct solution...
	mxAlgebra(cbind(-1.98*X1 + 0.2, 200*X2), name="grad", dimnames=list(NULL,c("x1","x2"))),
	mxFitFunctionAlgebra(algebra="BukinN2Func")#,gradient="grad")
)
set.seed(241216); m5 <- mxTryHard(model=m5,greenOK=T,checkHess=T)
omxCheckTrue(m5$output$status$code %in% c(0,1))
# Should still be TRUE if the sole free parameter has an active bound:
omxCheckTrue(m5$output$infoDefinite)

mxOption(reset=TRUE)
