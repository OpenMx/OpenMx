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

#options(error = browser)
require(OpenMx)
data <- mxData(type = 'raw', matrix(".", 3, 3, dimnames = list(NULL,c('a','b','c'))))
covariance <- mxMatrix('Symm', 3, 3, values = c(1:6), name = 'cov')
means <- mxMatrix('Full', 1, 3, values = c(1:3), name = 'means')
objective <- mxExpectationNormal('cov', 'means')
model <- mxModel('model', objective, covariance, means, data, mxFitFunctionML())
omxCheckError(mxRun(model), paste("The data object", omxQuotes("model.data"),
	"contains an observed matrix that is not of type 'double'"))
	
# Define a model
model <- mxModel()
model <- mxModel(model, mxMatrix("Full", values = c(0,0.2,0,0), name="A", nrow=2, ncol=2))
model <- mxModel(model, mxMatrix("Symm", values = c(0.8,0,0,0.8), name="S", nrow=2, ncol=2, free=TRUE))
model <- mxModel(model, mxMatrix("Iden", name="F", nrow=2, ncol=2, dimnames = list(c('a','b'), c('a','b'))))

model[["A"]]$free[2,1] <- TRUE
model[["S"]]$free[2,1] <- FALSE
model[["S"]]$free[1,2] <- FALSE
model[["S"]]$labels[1,1] <- "apple"
model[["S"]]$labels[2,2] <- "banana"

# Bounds must be added after all the free parameters are specified
model <- mxModel(model, mxBounds(c("apple", "banana"), 0.001, NA))

# Define the objective function
objective <- mxExpectationRAM("A", "S", "F")

# Define the observed covariance matrix
covMatrix <- matrix( c(0.77642931, 0.39590663, 0.39590663, 0.49115615), 
	nrow = 2, ncol = 2, byrow = TRUE, dimnames = list(c('a','b'), c('a','b')))

data <- mxData(covMatrix, 'cov', numObs = 100)
data$numObs <- 100L

# Add the objective function and the data to the model
model <- mxModel(model, objective, data, mxFitFunctionML())

fit <- mxRun(model)

primaryKey <- c(1:4, 3L)
m1 <- mxModel("uniqueModel", type="RAM",
              latentVars = "ign",
              mxData(type="raw", observed=data.frame(key=primaryKey), primaryKey = "key"),
              mxPath("one", "ign"),
              mxPath("ign", arrows=2))
omxCheckError(mxRun(m1), "uniqueModel.data: primary keys are not unique (examine rows with key=3)")

bad <- mxModel("bad", type="RAM",
              latentVars = "ign",
	       mxPath("one", "ign"),
	       mxPath("ign", arrows=2),
	       mxData(data.frame(key=1), 'raw', primaryKey="key"))
omxCheckError(mxRun(bad), "bad.data: primary key must be an integer or factor column in raw observed data")
