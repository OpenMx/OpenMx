#
#   Copyright 2007-2010 The OpenMx Project
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


# -----------------------------------------------------------------------
# Program:  SimpleCheckpoint.R
#  Author: Timothy R. Brick
#    Date: 05 12 2010 
#
# OpenMx one factor matrix model demo for front page of website
# 
# Revision History
#   Steven M. Boker -- 07 30 2009 created as OneFactorMatrixDemo.R
#   Hermine Maes -- 02 22 2010 updated & reformatted
#   Timothy R. Brick -- 05 12 2010 added checkpointing
# -----------------------------------------------------------------------

require(OpenMx)

data(demoOneFactor)
manifestVars <- names(demoOneFactor)

factorModel <- mxModel("One Factor",
    mxMatrix(type="Full", nrow=5, ncol=1, values=0.2, free=TRUE, name="A", labels=letters[1:5]),
    mxMatrix(type="Symm", nrow=1, ncol=1, values=1, free=FALSE, name="L"),
    mxMatrix(type="Diag", nrow=5, ncol=5, values=1, free=TRUE, name="U"),
    mxAlgebra(expression=A %*% L %*% t(A) + U, name="R"),
    mxMLObjective(covariance="R", dimnames=manifestVars),
    mxData(observed=cov(demoOneFactor), type="cov", numObs=500)
)

directory <- tempdir()

factorModel <- mxOption(factorModel, "Checkpoint Directory", directory)
factorModel <- mxOption(factorModel, "Checkpoint Units", "iterations")
factorModel <- mxOption(factorModel, "Checkpoint Count", 10)

factorFit <- mxRun(factorModel, checkpoint = TRUE)

factorRestore <- mxRestore(factorModel, chkpt.directory = directory)

omxCheckCloseEnough(mxEval(A, factorFit), mxEval(A, factorRestore), 0.001)
omxCheckCloseEnough(mxEval(U, factorFit), mxEval(U, factorRestore), 0.001)

