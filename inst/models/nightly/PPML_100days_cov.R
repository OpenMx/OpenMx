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


# One hundred days model
# Each manifest represents a single variable as measured on a different day, days 1-100
# Three predicting latents. Over the trial, one represents a constant part,
# one represents a linear part, and one an exponential part

require(OpenMx)


# Random covariance matrix

manifests <- unlist(lapply(1:100, function(n) { paste("Day",n, sep="") } ))
latents <- c('Const', 'Lin', 'Exp')

amtT <- 100
amtStep <- amtT / 99
ConstLoadings <- rep(1,100)
LinLoadings <- unlist(lapply(1:100, function (n) { n*amtStep/(amtT-1) } ))
EXPFACTOR <- 0.22
ExpLoadings <- unlist(lapply(1:100, function (n) { -exp(-EXPFACTOR*n*amtStep) } ))

lambda <- cbind(ConstLoadings, LinLoadings, ExpLoadings)
dataLatents <- matrix(rnorm(3*3), 3, 3)
dataTest <- lambda %*% dataLatents %*% t(dataLatents) %*% t(lambda) + diag(rep(1,100))
dataTest <- (dataTest + t(dataTest)) / 2

dataRaw <- mvtnorm::rmvnorm(n=5000, rep(0,100), dataTest)
colnames(dataRaw) <- manifests

colnames(dataTest) <- manifests
rownames(dataTest) <- manifests

factorModel <- mxModel("One Hundred Days Model",
      type="RAM",
	  # Vars
      manifestVars = manifests,
      latentVars = latents,
      # A
	  mxPath(from='Const',  to=manifests,value=ConstLoadings,	free=FALSE),
	  mxPath(from='Lin', 	to=manifests,value=LinLoadings,		free=FALSE),
	  mxPath(from='Exp',	to=manifests,value=ExpLoadings,		free=FALSE),
      # S
	  mxPath(from=manifests, arrows=2,value=rep(1.0, 100), labels=rep("Res", 100)),
	  #mxPath(from=latents, arrows=2,values=1.0),
	  mxPath(from="Const", to=latents, arrows=2,values=1.0),
	  mxPath(from="Lin", to=latents, arrows=2,values=1.0),
	  mxPath(from="Exp", to=latents, arrows=2,values=1.0),
	  
	  # Means
#	  mxPath(from="one", to=latents, values=0, free=TRUE),
	  
	  # Data
      mxData(dataTest, type="cov", numObs=100)
	)

factorModelOut <- mxRun(imxPPML(factorModel, TRUE))
