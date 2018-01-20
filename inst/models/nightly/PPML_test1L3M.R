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


require(OpenMx)


# Seed
set.seed(4)

# Data covariance matrix
dataCov <- rbind(c(4,2,2),
				 c(2,4,2),
				 c(2,2,4))
dataMean <- c(2,-2,1)

# Variable names
manifests <- c('X','Y','Z')
latents <- c('G')

# Dimnames
colnames(dataCov) <- manifests
rownames(dataCov) <- manifests
names(dataMean) <- manifests

factorModel <- mxModel("Factor Model 1L3M",
      type="RAM",
      manifestVars = manifests,
      latentVars = latents,
	  
	  # Regression Loadings
      mxPath(from='G', to=manifests,value=c(1,2,3),free=FALSE),
	  
	  # Variances
      mxPath(from=manifests, arrows=2,value=c(1,1,1), labels=c('E1','E1','E1')), # error variance
      mxPath(from=latents, arrows=2,values=1.0, labels=c("VG")), # variances of the latents
	 
	  # Latent means vector	 
      mxPath(from="one", to=latents, arrows=1, values=0, labels=c("MG"), free=TRUE),

	  mxData(dataCov, mean=dataMean, type="cov",numObs=100)
	  )
# Should be able to test fake latents, but there's a bug blocking the foldout function from working on this model
imxPPML.Test.Battery(factorModel, testFakeLatents=FALSE, tolerances=c(NA, .0001, .0001)) # NA -> Don't check covariance data w/ means

