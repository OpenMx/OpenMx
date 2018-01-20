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

data(demoOneFactor)
factorModel <- mxModel(
    "One Factor",
    mxMatrix("Full", 5, 1, values=0.2,
	     free=TRUE, name="A"),
    mxMatrix("Symm", 1, 1, values=1,
	     free=FALSE, name="L"),
    mxMatrix("Diag", 5, 5, values=1,
	     free=TRUE, name="U"),
    mxAlgebra(A %*% L %*% t(A) + U, name="R"),
    mxExpectationNormal("R", dimnames = names(demoOneFactor)),
    mxFitFunctionML(),
    mxData(cov(demoOneFactor), type="cov", numObs=500),
    mxComputeSequence(list(
      mxComputeNumericDeriv(),
      mxComputeReportDeriv()))
    )
fitModel <- mxRun(factorModel)
fullH <- fitModel$output$hessian

omxCheckEquals(fitModel$compute$steps[[1]]$output$probeCount, 4 * (10^2+10))

limModel <- mxRun(mxModel(factorModel,
                          mxComputeSequence(list(
                            mxComputeNumericDeriv(verbose=0,
                                                  knownHessian=fullH[3:10,3:10]),
                            mxComputeReportDeriv()))))
omxCheckCloseEnough(limModel$output$hessian[,1:2], fullH[,1:2], 1e-3)

omxCheckEquals(limModel$compute$steps[[1]]$output$probeCount, 152)
