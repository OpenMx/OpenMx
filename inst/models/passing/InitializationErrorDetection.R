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

data <- matrix(-1, dimnames=list(c("x"), c("x")))
omxCheckError(mxData(data, type="cov", numObs = 50),
	paste("The observed covariance matrix is not",
	"a positive-definite matrix:\n",
	 "1 or more elements of eigen(covMatrix)$values  <= 0"))
data <- matrix(-1, dimnames=list(c("x"), c("x")))
omxCheckError(mxData(data, type="cor", numObs = 50),
	paste("The observed correlation matrix is not",
	"a positive-definite matrix"))

data <- matrix(1, dimnames=list(c("x"), c("x")))
model <- mxModel('ErrorModel', 
    mxMatrix("Full", 1, 1, F, 1, name="cov", 
                dimnames=list(c("x"), c("x"))),
    mxData(data, type="cov", numObs = 50),
    mxFitFunctionML(),mxExpectationNormal("cov")
)
model$data$observed[1,1] <- -1
omxCheckError(mxRun(model), 
		"Observed Covariance Matrix is non-positive-definite.")
