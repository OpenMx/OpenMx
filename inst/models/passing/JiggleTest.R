#
#   Copyright 2007-2019 by the individuals mentioned in the source code history
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
manifests <- names(demoOneFactor)
latents <- c("G")
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
iniPars <- coef(factorModel)

#The largest start value is 1.0, so with scale=0.05 and dsn="runif",
#no perturbed parameter value can differ
#from its start value by more than 0.1:
pars2 <- imxJiggle(params=iniPars,lbounds=NA,ubounds=NA,dsn="runif",loc=1,scale=0.05)
omxCheckCloseEnough(iniPars,pars2,0.1)
mod2 <- mxJiggle(model=factorModel,scale=0.05)
omxCheckCloseEnough(iniPars,coef(mod2),0.1)

#Test classic:
mod3 <- mxJiggle(model=factorModel,classic=T)
omxCheckEquals(iniPars+0.1*(iniPars+0.5), coef(mod3))
