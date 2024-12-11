#
#   Copyright 2007-2020 by the individuals mentioned in the source code history
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


# http://openmx.ssri.psu.edu/issue/2014/05/memory-leak-when-running-ram-model-constraint

library(OpenMx)

#mxOption(NULL, "Default optimizer", 'SLSQP')

data(demoOneFactor)
manifests <- names(demoOneFactor)
latents <- c("G")
factorModelPath <- mxModel(
	"OneFactorPath",
	type="RAM",
	manifestVars = manifests,
	latentVars = latents,
	mxPath(from=latents, to=manifests,
				 labels=paste("l",1:5,sep="")),
	mxPath(from=manifests, arrows=2),
	mxPath(from=latents, arrows=2,
				 free=FALSE, values=1.0),
	mxData(cov(demoOneFactor), type="cov",
				 numObs=500),
	mxAlgebra(S[6,6],name="GV"),
	mxConstraint(GV-1==0,name="pointless"),
	mxConstraint(GV>0,name="morePointless"))
#factorModelPath <- mxOption(factorModelPath,"Checkpoint Directory","C:/Work/OpenMx_dev/")
#factorModelPath <- mxOption(factorModelPath,"Checkpoint Units","evaluations")
#factorModelPath <- mxOption(factorModelPath,"Checkpoint Count",1)
factorModelPath <- mxOption(factorModelPath,"Calculate Hessian","No")
factorModelPath <- mxOption(factorModelPath,"Standard Errors","No")
factorFit <- mxRun(factorModelPath, silent = TRUE)

omxCheckEquals(factorFit$output$status$code, 0)

omxCheckTrue(all(factorFit$output$constraintJacobian == 0))
omxCheckTrue(!("pointless" %in% names(factorFit$output$gradient)))
omxCheckTrue(!("morePointless" %in% names(factorFit$output$gradient)))

omxCheckCloseEnough(sqrt(sum(factorFit$output$gradient^2)), 0, .08)

