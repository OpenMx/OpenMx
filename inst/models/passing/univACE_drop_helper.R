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


# SCRIPT: univACE_drop_helper.R
# Timothy  Bates tim.bates$ed.ac.uk
# History:  Thu Oct  8 14:35:02 BST 2009
# OpenMx: http://www.openmx.virginia.com
##########################################

require(OpenMx)

# Prepare Data
data("twinData", package="OpenMx")
selVars <- c('bmi1','bmi2')
aceVars <- c("A1","C1","E1","A2","C2","E2")
mzfData <- subset(twinData, zyg==1, selVars)
dzfData <- subset(twinData, zyg==3, selVars)
cov(mzfData, use="pairwise.complete.obs")
cov(dzfData, use="pairwise.complete.obs")

# Fit ACE Model with RawData and Matrices Input
selVars <- c('x','y')
dataMZ <- matrix(c(1,.8,.8,1), nrow = 2, ncol=2, dimnames = list(selVars,selVars))
dataDZ <- matrix(c(1,.5,.5,1), nrow = 2, ncol=2, dimnames = list(selVars,selVars))

modelShare = mxModel("share", 
  mxMatrix("Full",values= .6, free=TRUE, labels='a1', nrow=1, ncol=1, name="a"),
  mxMatrix("Full",values= .6, free=TRUE, labels='c1', nrow=1, ncol=1, name="c"),
  mxMatrix("Full",values= .6, free=TRUE, labels='e1', nrow=1, ncol=1, name="e"),
  mxAlgebra(a * t(a), name="A"),
  mxAlgebra(c * t(c), name="C"),
  mxAlgebra(e * t(e), name="E"),
  mxAlgebra(name="MZcov",
            rbind(cbind(A+C+E, A+C),
                  cbind(A+C,   A+C+E)) 
  ),
  mxAlgebra(name="DZcov",
      rbind(cbind(A+C+E,     .5 %x%A+C),
            cbind(.5 %x%A+C,  A+C+E)) 
  )
)

modelMZ <- mxModel("modelMZ", 
  mxData(observed=dataMZ, type="cov", numObs = 100),
  mxFitFunctionML(),mxExpectationNormal(covariance = "share.MZcov",dimnames=selVars)  
)

modelDZ <- mxModel("modelDZ",
  mxData(observed=dataDZ, type="cov", numObs = 100),
  mxFitFunctionML(),mxExpectationNormal(covariance = "share.DZcov",dimnames=selVars)
)

model = mxModel("ACE", 
  modelShare, modelMZ, modelDZ,
  mxAlgebra(modelMZ.objective + modelDZ.objective, name="twin"),
  mxFitFunctionAlgebra("twin")
)
fit <- mxRun(model)
summary(fit)

# Now lets drop C using the helper
modelShare = omxSetParameters(modelShare, labels=c("c1"), free = FALSE, value = 0)

# rebuild and run the model
model = mxModel("ACE", 
  modelShare, modelMZ, modelDZ,
  mxAlgebra(modelMZ.objective + modelDZ.objective, name="twin"),
  mxFitFunctionAlgebra("twin")
)
fit <- mxRun(model)
summary(fit)
# hey presto: no c in the fit
