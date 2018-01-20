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


# -----------------------------------------------------------------------------
# Program:  SimpleCheckpoint.R
# Author: Timothy R. Brick
# Date: 2010.05.12 
#
# ModelType: Factor
# DataType: Continuous
# Field: None
#
# Purpose: 
#      OpenMx one factor matrix model demo for front page of website
# 
# RevisionHistory:
#      Steven M. Boker -- 2009.07.30 created as OneFactorMatrixDemo.R
#      Hermine Maes -- 2010.02.22 updated & reformatted
#      Timothy R. Brick -- 2010.05.12 added checkpointing
#      Ross Gore -- 2011.06.07	added Model, Data & Field metadata
# -----------------------------------------------------------------------------

require(OpenMx)
# Load Library
# -----------------------------------------------------------------------------

data(demoOneFactor)
# Prepare Data
# -----------------------------------------------------------------------------

manifestVars <- names(demoOneFactor)
# Prepare the manifest variables
# -----------------------------------------------------------------------------

factorModel <- mxModel("One Factor",
    mxMatrix(type="Full", nrow=5, ncol=1, values=0.2, free=TRUE, name="A", labels=letters[1:5]),
    mxMatrix(type="Symm", nrow=1, ncol=1, values=1, free=FALSE, name="L"),
    mxMatrix(type="Diag", nrow=5, ncol=5, values=1, free=TRUE, name="U"),
    mxAlgebra(expression=A %*% L %*% t(A) + U, name="R"),
    mxExpectationNormal(covariance="R", dimnames=manifestVars),
    mxFitFunctionML(),
    mxData(observed=cov(demoOneFactor), type="cov", numObs=500)
)
# Create an MxModel object using a matrix model specification
# -----------------------------------------------------------------------------

directory <- tempdir()
# Get a temporary directory for checkpointing
# -----------------------------------------------------------------------------

factorModel <- mxOption(factorModel, "Checkpoint Directory", directory)
factorModel <- mxOption(factorModel, "Checkpoint Units", "iterations")
factorModel <- mxOption(factorModel, "Checkpoint Count", 10)
# Prepare the model for the checkpointing
# -----------------------------------------------------------------------------

factorFit <- mxRun(factorModel, checkpoint = TRUE)
# Fit the model to the observed covariances
# -----------------------------------------------------------------------------

factorRestore <- mxRestore(factorModel, chkpt.directory = directory)
# Load the last saved state from the checkpoint file 
# -----------------------------------------------------------------------------

omxCheckCloseEnough(mxEval(A, factorFit), mxEval(A, factorRestore), 0.001)
omxCheckCloseEnough(mxEval(U, factorFit), mxEval(U, factorRestore), 0.001)
# Compare non-checkpointed results to checkpointed results 
# -----------------------------------------------------------------------------
