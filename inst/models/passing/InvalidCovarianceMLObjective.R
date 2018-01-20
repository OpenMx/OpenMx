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
cov <- mxMatrix('Full', 2, 2, values = c(0,1,1,0), name = 'cov', free = c(FALSE,TRUE,TRUE,FALSE))
objective <- mxExpectationNormal('cov', dimnames = c('a','b'))
identity <- diag(2)
dimnames(identity) <- list(c('a','b'),c('a','b'))
data <- mxData(identity, 'cov', numObs = 10)
model <- mxModel('model', cov, objective, data,  mxFitFunctionML())
ign <- omxCheckWarning(omxCheckError(mxRun(model), "The job for model 'model' exited abnormally with the error message: fit is not finite (model.fitfunction: stan::prob::multi_normal_sufficient_log: LDLT_Factor of covariance parameter is not positive definite.  last conditional variance is 0.)"),
		       "In model 'model' Optimizer returned a non-zero status code 10. Starting values are not feasible. Consider mxTryHard()")

dimnames(identity) <- list(c('a','c'),c('a','c'))
data <- mxData(identity, 'cov', numObs = 10)
model <- mxModel('model', cov, objective, data,  mxFitFunctionML())
ign <- omxCheckError(mxRun(model), "The dimnames for the expected covariance matrix ('a' and 'b') and the observed covariance matrix ('a' and 'c') in the Normal expectation function in model 'model' are not identical.")
