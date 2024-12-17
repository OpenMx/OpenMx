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

# The OpenMx backend should be able to recognize when active box constraints
# are making the Hessian spuriously appear non-PD, and NOT return status
# code 5 in such cases.

require(OpenMx)

data(demoOneFactor)

# Front-page model, nothing fancy:
factorModel <- mxModel(
	"One Factor",
	mxMatrix("Full", 5, 1, values=0.8, free=TRUE, name="A"),
	mxMatrix("Symm", 1, 1, values=1, free=FALSE, name="L"),
	mxMatrix("Diag", 5, 5, values=1, free=TRUE, name="U"),
	mxAlgebra(A %*% L %*% t(A) + U, name="R"),
	mxExpectationNormal(covariance = "R", dimnames = names(demoOneFactor)),
	mxFitFunctionML(),
	mxData(cov(demoOneFactor), type="cov", numObs=500)
)
summary(factorModelFit <- mxRun(factorModel))
omxCheckEquals(factorModelFit$output$status$code,0)
omxCheckTrue(factorModelFit$output$infoDefinite)

if (mxOption(NULL, "Default optimizer") == "SLSQP") { #<--TODO
	
	factorModel2 <- mxModel(
		"One Factor",
		# Notice the lbound:
		mxMatrix("Full", 5, 1, values=0.8, lbound=c(0.8,rep(NA,4)), free=TRUE, name="A"), 
		mxMatrix("Symm", 1, 1, values=1, free=FALSE, name="L"),
		mxMatrix("Diag", 5, 5, values=1, free=TRUE, name="U"),
		mxAlgebra(A %*% L %*% t(A) + U, name="R"),
		mxExpectationNormal(covariance = "R", dimnames = names(demoOneFactor)),
		mxFitFunctionML(),
		mxData(cov(demoOneFactor), type="cov", numObs=500)
	)
	fmf2 <- mxRun(factorModel2)
	summary(fmf2)
	omxCheckEquals(fmf2$output$status$code,0)
	omxCheckTrue(fmf2$output$infoDefinite)
	
	#Custom compute plan:
	plan <- omxDefaultComputePlan()
	plan$steps <- list(GD=plan$steps$GD,RE=plan$steps$RE)
	factorModel3 <- mxModel(
		"One Factor",
		plan,
		mxMatrix("Full", 5, 1, values=0.8, lbound=c(0.8,rep(NA,4)), free=TRUE, name="A"),
		mxMatrix("Symm", 1, 1, values=1, free=FALSE, name="L"),
		mxMatrix("Diag", 5, 5, values=1, free=TRUE, name="U"),
		mxAlgebra(A %*% L %*% t(A) + U, name="R"),
		mxExpectationNormal(covariance = "R", dimnames = names(demoOneFactor)),
		mxFitFunctionML(),
		mxData(cov(demoOneFactor), type="cov", numObs=500)
	)
	fmf3 <- mxRun(factorModel3)
	summary(fmf3)
	omxCheckEquals(fmf3$output$status$code,0)
	omxCheckTrue(!length(fmf3$output$infoDefinite))
	#^^^Because there was no MxComputeHessianQuality step in the compute plan.
	
	factorModel4 <- mxModel(
		"One Factor",
		# Note first loading fixed to 0.8:
		mxMatrix("Full", 5, 1, values=0.8, free=c(F,rep(T,4)), name="A"),
		mxMatrix("Symm", 1, 1, values=1, free=FALSE, name="L"),
		mxMatrix("Diag", 5, 5, values=1, free=TRUE, name="U"),
		mxAlgebra(A %*% L %*% t(A) + U, name="R"),
		mxExpectationNormal(covariance = "R", dimnames = names(demoOneFactor)),
		mxFitFunctionML(),
		mxData(cov(demoOneFactor), type="cov", numObs=500)
	)
	fmf4 <- mxRun(factorModel4)
	summary(fmf4)
	fmf4$output$standardErrors
	#^^^The standard errors should NOT be considered valid. 
	#It's "cheating" to fix an unknown parameter that, when freed, has an active bound at the MLE.
	
	fmf2.a <- mxTryHard(model=factorModel2,greenOK=T,checkHess=T)
	summary(fmf2.a)
	omxCheckEquals(fmf2.a$output$status$code,0)
	omxCheckTrue(fmf2.a$output$infoDefinite)
	
	fmf3.a <- mxTryHard(model=factorModel3,greenOK=T,checkHess=T)
	omxCheckEquals(fmf3.a$output$status$code,0)
	omxCheckTrue(!length(fmf3.a$output$infoDefinite))
	
}
