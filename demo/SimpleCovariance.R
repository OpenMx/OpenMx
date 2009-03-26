#
#   Copyright 2007-2009 The OpenMx Project
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
model <- mxModel(model, mxMatrix("Full", c(0,0.2,0,0), name="A", nrow=2, ncol=2))
model <- mxModel(model, mxMatrix("Full", c(0.8,0,0,0.8), name="S", nrow=2, ncol=2, free=TRUE))
model <- mxModel(model, mxMatrix("Full", c(1,0,0,1), name="F", nrow=2, ncol=2))
model <- mxModel(model, mxBounds(c("apple", "banana"), 0, NA))

model[["A"]]@specification[2,1] <- NA
model[["S"]]@specification[2,1] <- 0
model[["S"]]@specification[1,2] <- 0
model[["S"]]@specification[1,1] <- "apple"
model[["S"]]@specification[2,2] <- "banana"

# Define the objective function
objective <- mxRAMObjective("A", "S", "F", "objective")

# Define the observed covariance matrix
covMatrix <- matrix( c(0.77642931, 0.39590663, 0.39590663, 0.49115615), nrow = 2, ncol = 2, byrow = TRUE)

# Add the objective function and the data to the model
model <- mxModel(model, objective, covMatrix)

# Run the job
model <- mxJobRun(model)

print(model@output)
