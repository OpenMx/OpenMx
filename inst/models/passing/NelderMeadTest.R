#
#   Copyright 2007-2017 The OpenMx Project
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
foo <- mxComputeNelderMead()
plan <- omxDefaultComputePlan()
plan$steps <- list(foo,plan$steps$RE)

data(demoOneFactor)
manifests <- names(demoOneFactor)
latents <- c("G")
factorModel <- mxModel("OneFactor",
											 type="RAM",
											 plan,
											 manifestVars = manifests,
											 latentVars = latents,
											 mxPath(from=latents, to=manifests),
											 mxPath(from=manifests, arrows=2),
											 mxPath(from=latents, arrows=2,
											 			 free=FALSE, values=1.0),
											 mxData(cov(demoOneFactor), type="cov",
											 			 numObs=500))
omxCheckError(mxRun(factorModel),"NelderMeadOptimizerContext::invokeNelderMead() : so far, so good")
