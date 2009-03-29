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

# This alternative version of SimpleCovariance spreads out the definitions
#   of the mx objects as much as possible.

############################################################
# Define the observed covariance matrix

covMatrix <- matrix( c(0.77642931, 0.39590663, 
                       0.39590663, 0.49115615), 
                    nrow = 2, ncol = 2, byrow = TRUE)
theData <- mxData(covMatrix, 'cov', numObs = 100, name = 'covariance matrix')

############################################################
# Define the model matrices

# First, define the asymmetric paths and starting values, the A matrix

AVal  <- matrix(c(0.0, 0.0,
                  0.2, 0.0), 
                nrow=2, ncol=2, byrow=T)
ASpec <- matrix(c(  0,   0,
                   NA,   0), 
                nrow=2, ncol=2, byrow=T)
theAMatrix <- mxMatrix(values=AVal, specification=ASpec, name="A")

# Next, define the symmetric paths and starting values, the S matrix

SVal  <- matrix(c(0.8, 0.0,
                  0.0, 0.8), 
                nrow=2, ncol=2, byrow=T)
SSpec <- matrix(c("apple",   0,
                    0,   "banana"), 
                nrow=2, ncol=2, byrow=T)
theSMatrix <- mxMatrix(values=SVal, specification=SSpec, name="S")

# Finally, define the filter matrix as a 2x2 identity matrix

theFMatrix <- mxMatrix("Iden", name="F", nrow=2, ncol=2)

############################################################
# Define the bounds

theBounds1 <- mxBounds(c("apple", "banana"), 0.001, NA, name="theBounds1")

############################################################
# Define the model

modelAlt <- mxModel(theAMatrix, theSMatrix, theFMatrix, theBounds1,
                    mxRAMObjective(name="ramObjective"), theData)

modelAlt[["A"]]
modelAlt[["S"]]
modelAlt[["F"]]
modelAlt[["theBounds1"]]
modelAlt@data
modelAlt@objective

############################################################
# Run the job

modelAltOut <- mxJobRun(modelAlt)

print(modelAltOut@output)

modelAltOut[["A"]]
modelAltOut[["S"]]
modelAltOut[["F"]]
modelAltOut@data
modelAltOut@objective


