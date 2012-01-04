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
# Program: OneFactorPathDemo.R  
# Author: Steve Boker
# Date: 2009.08.01 
#
# ModelType: Factor
# DataType: Continuous
# Field: None
#
# Purpose: 
#      OpenMx one factor path model demo for front page of website
# 
# RevisionHistory:
#      Hermine Maes -- 2010.02.22 updated & reformatted
#      Ross Gore -- 2011.06.06	added Model, Data & Field metadata
# -----------------------------------------------------------------------------

require(OpenMx)
# Load Library
# -----------------------------------------------------------------------------

data(demoOneFactor)
# Prepare Data
# -----------------------------------------------------------------------------

manifests <- names(demoOneFactor)
# Prepare Manifests Data
# -----------------------------------------------------------------------------

latents <- c("G")
# Prepare Latents Data
# -----------------------------------------------------------------------------

factorModel <- mxModel("One Factor", 
    type="RAM",
    manifestVars=manifests, 
    latentVars=latents,
    mxPath(from=latents, to=manifests),
    mxPath(from=manifests, arrows=2),
    mxPath(from=latents, arrows=2, free=FALSE, values=1.0),
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