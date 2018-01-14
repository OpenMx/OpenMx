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

if (mxOption(NULL, "Default optimizer") == "SLSQP") {
  mxOption(NULL, "Optimality tolerance", "6e-10")
}

set.seed(654684)
data(demoOneFactor)

manifests <- names(demoOneFactor)
latents <- c("factor")
nVar <- length(manifests)
nFac <- length(latents)

factorModel <- mxModel("One Factor ML",
    mxData(cov(demoOneFactor), type="cov", means=colMeans(demoOneFactor), numObs=500),
    mxMatrix("Full", 1, nVar, free=T, values=colMeans(demoOneFactor), name="M"),
    mxMatrix("Full", nVar, nFac, free=T, values=.2, ubound=1, name="A"),
    mxMatrix("Diag", nVar, nVar, free=T, values=1, lbound=.0001, ubound=10, name="D"),
    mxAlgebra(A%*%t(A) + D, name="C"),
    mxAlgebra(sqrt(diag2vec(C)),name="P"),
    mxFitFunctionML(),mxExpectationNormal("C", "M", dimnames=manifests),
    mxCI(c("P"))
)
factorFitCI <- mxTryHard(factorModel, fit2beat=-2863, extraTries=3, intervals=TRUE)
factorSummCI <- summary(factorFitCI)

#print(factorFitCI$output$computes[[2]])

( ci <- factorFitCI$output$confidenceIntervals )
#cat(deparse(round(ci[,'lbound'], 4)))
omxCheckCloseEnough(ci[,'estimate'], c(0.4455, 0.54, 0.6115, 0.7302, 0.8186), 1e-3)
omxCheckCloseEnough(ci[,'lbound'] - c(0.4174, 0.5082, 0.5755, 0.6872, 0.7704), rep(0,5), 4e-2)
if (mxOption(NULL, "Default optimizer") != "NPSOL") {
  omxCheckCloseEnough(ci[,'ubound'], c(0.4747, 0.5754, 0.6516, 0.7781, 0.8723), 4e-2)
}

factorModelRaw <- mxModel(factorFitCI,
                          mxData(demoOneFactor, type="raw"),
                          mxFitFunctionML(),
                          mxExpectationNormal("C", "M", dimnames=manifests),
                          name = "One Factor FIML")
factorFitRawCI <- mxRun(factorModelRaw, intervals=TRUE, suppressWarnings = TRUE)
factorSummRawCI <- summary(factorFitRawCI)

ci <- factorFitRawCI$output$confidenceIntervals
omxCheckCloseEnough(ci[,'estimate'], c(0.4455, 0.54, 0.6115, 0.7302, 0.8186), 1e-3)
omxCheckCloseEnough(ci[,'lbound'] - c(0.4193, 0.5082, 0.5755, 0.6872, 0.7704), rep(0,5), 1e-3)
omxCheckCloseEnough(ci[,'ubound'], c(0.4747, 0.5754, 0.6516, 0.7781, 0.8723), 1e-3)

# Compare to original MX Estimates
#          5  Confidence intervals requested in group            1
# Matrix Element Int.      Estimate         Lower         Upper  Lfail Ufail
# P   1   1   1  95.0         0.4465        0.4192        0.4749 0 1   0 1
# P   1   1   2  95.0         0.5412        0.5081        0.5756 0 1   6 1
# P   1   1   3  95.0         0.6131        0.5753        0.6518 0 1   4 1
# P   1   1   4  95.0         0.7322        0.6870        0.7783 6 1   0 0
# P   1   1   5  95.0         0.8209        0.7702        0.8726 0 1   0 0
