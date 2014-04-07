#
#   Copyright 2007-2014 The OpenMx Project
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


require(OpenMx)
data(demoOneFactor)
setwd(tempdir())
factorModel <- mxModel("One_Factor",
                       mxMatrix("Full", 5, 1, values=0.2,
                                free=TRUE, name="A"),
                       mxMatrix("Symm", 1, 1, values=1,
                                free=FALSE, name="L"),
                       mxMatrix("Diag", 5, 5, values=1,
                                free=TRUE, name="U"),
                       mxAlgebra(A %*% L %*% t(A) + U, name="R"),
                       #mxExpectationNormal(covariance="R",dimnames=names(demoOneFactor)),
                       #mxFitFunctionML(),
                       mxMLObjective("R", dimnames = names(demoOneFactor)),
                       mxData(cov(demoOneFactor), type="cov", numObs=500))
factorModel <- mxOption(factorModel,"Checkpoint Units",'iterations')
factorModel <- mxOption(factorModel,"Checkpoint Count",1)
factorFit <- mxRun(factorModel,checkpoint=T)
checkpointdat <- read.table("One_Factor.omx",header=T,as.is=T)
omxCheckTrue(!is.null(checkpointdat))
omxCheckTrue(nrow(checkpointdat)>1)

factorModel2 <- mxOption(factorModel,"Default optimizer","NPSOL")
factorFit2 <- mxRun(factorModel2,checkpoint=T)
checkpointdat <- read.table("One_Factor.omx",header=T,as.is=T)
omxCheckTrue(!is.null(checkpointdat))
omxCheckTrue(nrow(checkpointdat)>1)

