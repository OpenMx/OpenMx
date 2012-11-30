#
#   Copyright 2007-2012 The OpenMx Project
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
# Program: OneFactorMatrixDemo.R  
# Author: Steve Boker
# Date: 2009.07.30 
#
# ModelType: Factor
# DataType: Continuous
# Field: None
#
# Purpose:
#      OpenMx one factor matrix model demo for front page of website
# 
# RevisionHistory:
#      Hermine Maes -- 2010.02.22 updated & reformatted
#      Ross Gore -- 2011.06.06 added Model, Data & Field metadata
# -----------------------------------------------------------------------------

require(OpenMx)
# Load Library
# -----------------------------------------------------------------------------

data(demoOneFactor)
# Prepare Data
# -----------------------------------------------------------------------------


manifestVars <- names(demoOneFactor)
# Prepare Manifests Data
# -----------------------------------------------------------------------------

factorModel <- mxModel("One Factor",
    mxMatrix(type="Full", nrow=5, ncol=1, values=0.2, free=TRUE, name="A"),
    mxMatrix(type="Symm", nrow=1, ncol=1, values=1, free=FALSE, name="L"),
    mxMatrix(type="Diag", nrow=5, ncol=5, values=1, free=TRUE, name="U"),
    mxAlgebra(expression=A %*% L %*% t(A) + U, name="R"),
    mxFitFunctionML(),mxExpectationNormal(covariance="R", dimnames=manifestVars),
    mxData(observed=cov(demoOneFactor), type="cov", numObs=500)
)
# Create an MxModel object
# -----------------------------------------------------------------------------

factorFit <- mxRun(factorModel)
# Fit the model to the observed covariances with mxRun
# -----------------------------------------------------------------------------

summary(factorFit)
# Print a summary of the results
# -----------------------------------------------------------------------------
