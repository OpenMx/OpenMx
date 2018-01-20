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

varNames <- c('x')

data1 <- mxData(matrix(1, dimnames = list(varNames,varNames)), type="cov", numObs=100)
data2 <- mxData(matrix(2, dimnames = list(varNames,varNames)), type="cov", numObs=100)

mat1 <- mxMatrix("Full",2,free=TRUE, nrow=1, ncol=1, name="mat1")
mat2 <- mxMatrix("Full",1,free=TRUE, nrow=1, ncol=1, name="mat2")

obj1 <- mxExpectationNormal("mat1", dimnames = varNames)
obj2 <- mxExpectationNormal("mat2", dimnames = varNames)

model1 <- mxModel("model1", data1, mat1, obj1, mxFitFunctionML())
model2 <- mxModel("model2", data2, mat2, obj2, mxFitFunctionML())

output1 <- mxRun(model1, suppressWarnings = TRUE)
output2 <- mxRun(model2, suppressWarnings = TRUE)

output1$output
output2$output

alg <- mxAlgebra(model1.objective + model2.objective, name="alg")
obj <- mxFitFunctionAlgebra("alg")

model <- mxModel("both", alg, obj, model1, model2)
model <- mxRun(model, suppressWarnings = TRUE)

omxCheckCloseEnough(model$output$estimate, .99 * c(1, 2), 0.001)

# Also use the multigroup fit function
mft <- mxFitFunctionMultigroup(c("model1", "model2"))
mgmodel <- mxModel("mgboth", mft, model1, model2)
mgmodel <- mxRun(mgmodel, suppressWarnings = TRUE)

#check estimates
omxCheckCloseEnough(omxGetParameters(model), omxGetParameters(mgmodel), 1e-3)

#check standard errors
omxCheckCloseEnough(summary(model)$parameters[,6], summary(mgmodel)$parameters[,6], 1e-3)

mgmodel$output$hessian
mgmodel$output$calculatedHessian
mgmodel$output$hessianCholesky


omxCheckError(mxCheckIdentification(model), "Identification check is not possible for models with 'MxFitFunctionAlgebra', 'MxFitFunctionRow', and 'MxFitFunctionR' fit functions.
 If you have a multigroup model, use mxFitFunctionMultigroup.")

