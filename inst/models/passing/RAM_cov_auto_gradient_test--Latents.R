#
#   Copyright 2007-2023 by the individuals mentioned in the source code history
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
#mxOption(NULL,"Default optimizer","SLSQP")
#mxOption(NULL,"Verify level",-1)
data(demoOneFactor)
manifests <- names(demoOneFactor)
latents <- c("G")
mxOption(NULL,"Analytic gradients","Yes")
factorModel <- mxModel(
	"One Factor",
	type="RAM",
	manifestVars = manifests,
	latentVars = latents,
	mxPath(from=latents, to=manifests,values=0.8),
	mxPath(from=manifests, arrows=2,values=1),
	mxPath(from=latents, arrows=2,
				 free=FALSE, values=1.0),
	mxData(cov(demoOneFactor), type="cov",
				 numObs=500)
)
f1 <- mxRun(factorModel)
mxOption(NULL,"Analytic gradients","No")
f2 <- mxRun(factorModel)
omxCheckCloseEnough(coef(f1)-coef(f2),rep(0,10),2e-6)
omxCheckCloseEnough(f1$output$gradient-f2$output$gradient,rep(0,10),1.5e-2)
omxCheckCloseEnough(f1$output$fit-f2$output$fit,0,5e-8)

#Using analytic derivatives should be faster:
omxCheckTrue(f1$output$iterations <= f2$output$iterations)
f1$output$iterations; f2$output$iterations
omxCheckTrue(f1$output$evaluations < f2$output$evaluations)
f1$output$evaluations ; f2$output$evaluations
omxCheckTrue(summary(f1)$wallTime < summary(f2)$wallTime)
summary(f1)$wallTime ; summary(f2)$wallTime
