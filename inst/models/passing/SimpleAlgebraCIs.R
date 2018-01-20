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

# give SLSQP plenty of time to converge
mxOption(NULL, "major iterations", 3000)

data(demoOneFactor)

manifests <- names(demoOneFactor)
latents <- c("factor")
nManifest <- length(manifests)
nVars <- nManifest + length(latents)

factorModel <- mxModel("One Factor", type="RAM",
    manifestVars = manifests,
    latentVars = latents,
    mxPath(from=latents, to=manifests, free=c(FALSE,TRUE,TRUE,TRUE,TRUE), values=1),
    mxPath(from=manifests, arrows=2, lbound=.0001),
    mxPath(from=latents, arrows=2, free=TRUE, values=1.0),
    mxPath(from="one", to=manifests, arrows=1, free=T, values=colMeans(demoOneFactor)),
    mxData(demoOneFactor, type="raw"),
    mxMatrix("Iden", nrow=nVars, name="I"),
    mxMatrix("Full", free=FALSE, values=diag(nrow=nManifest, ncol=nVars), name="Eff"),
    mxAlgebra(Eff%*%solve(I-A), name="Z"),
    mxAlgebra(Z%*%S%*%t(Z), name="C"),
    mxAlgebra(sqrt(diag2vec(C)), name="P"),
    mxCI(c("P"))
)
factorFit <- mxRun(factorModel, intervals=FALSE)
omxCheckCloseEnough(factorFit$output$fit, 934.095, .01)

factorBoot <- mxBootstrap(factorFit, 100L, OK=c("OK","OK/green", "nonzero gradient/red"))
omxCheckError(mxBootstrapEval(P, factorFit, compute=T),
	      "Compute plan MxComputeSequence found in model 'One Factor' instead of MxComputeBootstrap. Have you run this model through mxBootstrap already?")

if (mxOption(NULL, 'Default optimizer') != "SLSQP") {ctype = 'none'} else {ctype = 'ineq'}

factorFitCI <- mxRun(mxModel(factorFit, mxComputeConfidenceInterval(plan=mxComputeGradientDescent(), constraintType = ctype)), suppressWarnings = TRUE)
factorSummCI <- summary(factorFitCI)
summary(factorFitCI)

set.seed(42L)
bci <- mxBootstrapEval(P, factorBoot, bq=c(.025,.975))
print(bci)
omxCheckCloseEnough(factorFitCI$output$confidenceIntervals[,'lbound'] - bci[,"2.5%"],
                    rep(0,5), .05)
omxCheckCloseEnough(factorFitCI$output$confidenceIntervals[,'ubound'] - bci[,"97.5%"],
                    rep(0,5), .03)

omxCheckCloseEnough(coef(factorFit), coef(factorFitCI))
omxCheckCloseEnough(factorFit$output$fit, factorFitCI$output$fit, 0)
omxCheckCloseEnough(mxEval(Z, factorFit), mxEval(Z, factorFitCI))
omxCheckCloseEnough(mxEval(C, factorFit), mxEval(C, factorFitCI))
omxCheckCloseEnough(mxEval(P, factorFit), mxEval(P, factorFitCI))
omxCheckCloseEnough(mxEval(objective, factorFit)[1,1], mxEval(objective, factorFitCI)[1,1])

if (0) {
  options(digits=12)
  print(factorFitCI$output$fit)
  print(factorFitCI$output$computes[[2]])
}

ci <- factorFitCI$output$confidenceIntervals
#print(ci)
#cat(deparse(round(ci[,'ubound'],4)))
omxCheckCloseEnough(ci[,'estimate'], c(0.4456, 0.5401, 0.6116, 0.7302, 0.8187), .001)

lbound <- c(0.406, 0.485, 0.553, 0.6872, 0.769)
ubound <- c(0.4747, 0.5754, 0.6516, 0.778, 0.8723)

if (mxOption(NULL, 'Default optimizer') == "CSOLNP") {
	omxCheckCloseEnough(ci[,'lbound'], lbound, .06)
} else {
	omxCheckCloseEnough(ci[,'lbound'], lbound, .03)
}

if (mxOption(NULL, 'Default optimizer') != "NPSOL") {
	# NPSOL needs to get slightly closer to the MLE to nail all of these
	omxCheckCloseEnough(ci[,'ubound'], ubound, .06)
}

factorParallel <- omxParallelCI(factorFit)
pci <- factorParallel$output$confidenceIntervals
omxCheckCloseEnough(pci[,'lbound'], lbound, .06)
omxCheckCloseEnough(pci[,'ubound'], ubound, .06)
