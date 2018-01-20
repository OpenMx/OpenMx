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

manifests <- names(demoOneFactor)
latents <- c("factor")
nVar <- length(manifests)
nFac <- length(latents)

factorModel <- mxModel("One Factor ML",
    mxData(cov(demoOneFactor), type="cov", numObs=500),
    mxMatrix("Full", 1, nVar, free=T, values=colMeans(demoOneFactor), name="M"),
    mxMatrix("Full", nVar, nFac, free=T, values=.2, name="A"),
    mxMatrix("Diag", nVar, nVar, free=T, values=1, lbound=.0001, name="D"),
    mxAlgebra(A %*% t(A) + D, name="C"),
    mxAlgebra(sqrt(diag2vec(C)),name="P"),
    mxFitFunctionML(),mxExpectationNormal("C", "M", dimnames=manifests))

omxCheckError(mxRun(factorModel), 
	paste("In model 'One Factor ML' the Normal expectation function",
		"contains an expected means vector but the model is missing",
		"some data for the observed means."))

factorModel <- mxModel("One Factor ML",
    mxData(cov(demoOneFactor), means = colMeans(demoOneFactor), type="cov", numObs=500),
    mxMatrix("Full", nVar, 1, free=T, values=colMeans(demoOneFactor), name="M"),
    mxMatrix("Full", nVar, nFac, free=T, values=.2, name="A"),
    mxMatrix("Diag", nVar, nVar, free=T, values=1, lbound=.0001, name="D"),
    mxAlgebra(A %*% t(A) + D, name="C"),
    mxAlgebra(sqrt(diag2vec(C)),name="P"),
    mxFitFunctionML(),mxExpectationNormal("C", "M", dimnames=manifests))

omxCheckError(mxRun(factorModel),
	paste("The expected means vector associated with",
		"the expectation function in model 'One Factor ML'",
		"is not a 1 x n matrix. It has dimensions 5 x 1."))

factorModel <- mxModel("One Factor FIML",
    mxData(demoOneFactor, type="raw"),
    mxMatrix("Full", nVar, 1, free=T, values=colMeans(demoOneFactor), name="M"),
    mxMatrix("Full", nVar, nFac, free=T, values=.2, name="A"),
    mxMatrix("Diag", nVar, nVar, free=T, values=1, lbound=.0001, name="D"),
    mxAlgebra(A %*% t(A) + D, name="C"),
    mxAlgebra(sqrt(diag2vec(C)),name="P"),
    mxFitFunctionML(),mxExpectationNormal("C", "M", dimnames=manifests))

omxCheckError(mxRun(factorModel),
	paste("The expected means vector associated with",
		"the expectation function in model 'One Factor FIML'",
		"is not a 1 x n matrix. It has dimensions 5 x 1."))
