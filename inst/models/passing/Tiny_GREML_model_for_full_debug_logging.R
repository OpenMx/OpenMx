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

set.seed(1234)
dat <- cbind(rnorm(5),rep(1,5))
colnames(dat) <- c("y","x")
dat[5,1] <- NA

ge <- mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=FALSE)
gff <- mxFitFunctionGREML(dV=c(ve="I"))
plan <- mxComputeSequence(steps=list(
	mxComputeNewtonRaphson(fitfunction="fitfunction"),
	mxComputeOnce('fitfunction', c('fit','gradient','hessian','ihessian')),
	mxComputeStandardError(),
	mxComputeReportDeriv(),
	mxComputeReportExpectation()
))

testmod <- mxModel(
	"GREMLtest",
	mxData(observed = dat, type="raw", sort=FALSE),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
	mxMatrix(type="Unit", nrow=5,ncol=1,name="Uni"),
	mxAlgebra(vec2diag(Uni), name="I"),
	mxAlgebra( vec2diag(Uni%x%Ve), name="V"),
	ge,
	gff,
	plan
)

testrun <- mxRun(testmod)

summary(testrun)
omxCheckEquals(mxEval(nrow(I),testrun,T), 5)
omxCheckEquals(mxEval(ncol(I),testrun,T), 5)

#mxGetExpected() knows how to drop rows and columns of the covariance matrix due to missing data; mxEvalByName() does not:
omxCheckEquals(nrow(mxEvalByName("V", testrun, compute=T)), 5)
omxCheckEquals(nrow(mxGetExpected(testrun, "covariance")), 4)

omxCheckTrue( all(abs(mxGetExpected(testrun,"means") - mean(dat[,1],na.rm=T)) < 1e-8) )
omxCheckTrue(is.na(mxGetExpected(testrun,"thresholds")))


testmod2 <- mxModel(
	"GREMLtest",
	mxData(observed = dat, type="raw", sort=FALSE),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
	mxMatrix(type="Unit", nrow=5,ncol=1,name="Uni"),
	mxAlgebra(vec2diag(Uni), name="I"),
	mxAlgebra( I%x%Ve, name="V"),
	ge,
	gff,
	plan
)

testrun2 <- mxRun(testmod2)

summary(testrun2)
omxCheckEquals(mxEval(nrow(I),testrun2,T), 5)
omxCheckEquals(mxEval(ncol(I),testrun2,T), 5)

