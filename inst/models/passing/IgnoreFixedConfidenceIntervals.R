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

# Add the objective function and the data to the model
model <- mxModel(model, objective, data, mxCI("F"), mxFitFunctionML())

print(model)

# Run the job
modelOut <- mxRun(model, intervals=TRUE)

expectedParameters <- c(0.5099, 0.7686, 0.2863)

omxCheckCloseEnough(expectedParameters, 
	modelOut$output$estimate, 
	epsilon = 10 ^ -4)

omxCheckEquals(length(summary(modelOut)$CI), 0) 

