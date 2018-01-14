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


# At one point we had a bug such that models could not be re-run without generating an error
# This tests that we don't get a regression on that bug

library(OpenMx)
data(myFADataRaw, package="OpenMx")
manifests = paste("x", 1:3, sep="")
latents   = c("G")
m1 <- mxModel("m1", type="RAM",
	manifestVars = manifests, latentVars   = latents,
	mxPath(from = latents, to = manifests),
	mxPath(from = manifests, arrows = 2), # manifest residuals 
	mxPath(from = latents, arrows = 2, free = FALSE, values = 1), # latent fixed at 1
	mxPath(from = "one", to = manifests, arrows = 1), # manifest means
	mxData(myFADataRaw[1:100, manifests], type = "raw")
)
diag(m1$S$lbound) <- .01
m1 = mxRun(m1)
m1 = mxRun(m1)
