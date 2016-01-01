#
#   Copyright 2007-2016 The OpenMx Project
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

model <- mxModel('model')
obsData <- matrix(c(1:10), 10, 1)
data <- mxData(obsData, "raw")
chol <- mxMatrix("Full", 2, 2, c(T,T,F,T), c(1,.2,0,1), name="Chol") 
expCov <- mxAlgebra(Chol %*% t(Chol), name="expCov") 
expMean <- mxMatrix("Full", 1, 2, T, c(0,0), "expMean")
objective <- mxExpectationNormal("expCov", "expMean")
foo <- mxMatrix('Full', 1, 1, name = 'foo')
bar <- mxMatrix('Full', 10, 1, name = 'bar')
model <- mxModel(model, data, chol, expCov, expMean, objective, foo, bar, mxFitFunctionML())
omxCheckIdentical(mxEval(objective + foo, model, TRUE), matrix(as.numeric(NA), 1, 1))
omxCheckError(mxEval(objective + bar, model, TRUE), "The following error occurred while evaluating the expression 'model.fitfunction + model.bar' in model 'model' : non-conformable arrays")
objective <- mxExpectationNormal("expCov", "expMean")
model <- mxModel(model, objective, mxFitFunctionML(vector = TRUE))
omxCheckError(mxEval(objective + foo, model, TRUE), "The following error occurred while evaluating the expression 'model.fitfunction + model.foo' in model 'model' : non-conformable arrays")
omxCheckIdentical(mxEval(objective + bar, model, TRUE), matrix(as.numeric(NA), 10, 1))
