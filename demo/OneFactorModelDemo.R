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
# Program: OneFactorMatrixDemo.R  
#  Author: Steve Boker
#    Date: 08 01 2009 
#
# OpenMx one factor matrix model demo from front page of website
# 
# Revision History
#   Hermine Maes -- 10 08 2009 updated & reformatted
# -----------------------------------------------------------------------

require(OpenMx)

data(demoOneFactor)

factorModel <- mxModel(name ="One Factor",
    mxMatrix(type="Full", nrow=5, ncol=1, free=T, values=0.2, name="A"),
    mxMatrix(type="Symm", nrow=1, ncol=1, free=T, values=1, name="L"),
    mxMatrix(type="Diag", nrow=5, ncol=5, free=T, values=1, name="U"),
    mxAlgebra(expression=A %*% L %*% t(A) + U, name="R"),
    mxMLObjective(covariance="R", dimnames=names(demoOneFactor)),
    mxData(observed=cov(demoOneFactor), type="cov", numObs=500))

factorModelFit <- mxRun(factorModel)
summary(factorModelFit)
