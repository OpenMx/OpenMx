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

#options(error = utils::recover)
require(OpenMx)

data(multiData1)

manifests <- c("x1", "y")

uniRegModelRaw <- mxModel("uniRegModelRaw",
    type="RAM",
    manifestVars=manifests,
    mxPath(from="x1", to="y", arrows=1, 
           free=TRUE, values=.2, labels="b1"),
    mxPath(from=manifests, 
           arrows=2, free=TRUE, values=.8, 
           labels=c("VarX1", "VarE")),
    mxPath(from="one", to=manifests, 
           arrows=1, free=TRUE, values=.1, 
           labels=c("MeanX1", "MeanY")),
    mxData(observed=multiData1, type="raw"),
    mxFitFunctionML(vector=TRUE)  # should fail
    )

varNames <- c('x')

data1 <- mxData(matrix(1, dimnames = list(varNames,varNames)), type="cov", numObs=100)
data2 <- mxData(matrix(2, dimnames = list(varNames,varNames)), type="cov", numObs=100)

mat1 <- mxMatrix("Full",2,free=TRUE, nrow=1, ncol=1, name="mat1")
mat2 <- mxMatrix("Full",1,free=TRUE, nrow=1, ncol=1, name="mat2")

obj1 <- mxExpectationNormal("mat1", dimnames = varNames)
obj2 <- mxExpectationNormal("mat2", dimnames = varNames)

model1 <- mxModel("model1", data1, mat1, obj1, mxFitFunctionML())
model2 <- mxModel("model2", data2, mat2, obj2, mxFitFunctionML())

#output1 <- mxRun(model1, suppressWarnings = TRUE)
#output2 <- mxRun(model2, suppressWarnings = TRUE)

#output1$output
#output2$output

omxCheckError(mxFitFunctionMultigroup(groups=c()),
              message="mxFitFunctionMultigroup: at least 1 fitfunction must be provided")

alg <- mxAlgebra(model1.objective + model2.objective, name="alg")
if (1) {
	obj <- mxFitFunctionMultigroup(paste("model", 1:2, sep=""))
	model <- mxModel("both", obj, model1, model2)
        model.est <- mxRun(model, suppressWarnings = TRUE)
        omxCheckCloseEnough(model.est$output$estimate, 99/100*c(1, 2), 0.001)
}
if (1) {
	obj <- mxFitFunctionMultigroup("both.alg")
	model <- mxModel("both", obj, model1, model2, alg)
        model.est <- mxRun(model, suppressWarnings = TRUE)
        omxCheckCloseEnough(model.est$output$estimate, 99/100*c(1, 2), 0.001)
}
if (1) {
	obj <- mxFitFunctionAlgebra("alg")
	model <- mxModel("both", obj, model1, model2, alg)
        model.est <- mxRun(model, suppressWarnings = TRUE)
        omxCheckCloseEnough(model.est$output$estimate, 99/100*c(1, 2), 0.001)
}
if (1) {
	obj <- mxFitFunctionMultigroup(c("uniRegModelRaw", paste("model", 1:2, sep="")))
	model <- mxModel(model="vector", obj, model1, model2, uniRegModelRaw)
	omxCheckError(mxRun(model, suppressWarnings = TRUE), "vector.fitfunction: cannot combine units Pr and -2lnL (from model1.fitfunction)")
}

model <- mxModel(model="fail", mxFitFunctionMultigroup("noExisto"))
omxCheckError(mxRun(model), "fail.fitfunction: cannot locate algebra/fitfunction 'noExisto'")
