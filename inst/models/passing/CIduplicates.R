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


library(OpenMx)
data(myFADataRaw, package="OpenMx")
manifests = paste0("x",1:3)
myFADataRaw = myFADataRaw[, manifests]
latents = c("G")
m1 <- mxModel("m1", type="RAM",
	manifestVars = manifests, latentVars   = latents,
	mxPath(from = latents, to = manifests),
	mxPath(from = manifests, arrows = 2, labels = paste0(manifests, "_resid")),
	mxPath(from = latents, arrows = 2, free = F, values = 1), # latents fixed at 1
	mxData(cov(myFADataRaw, use="complete"), type = "cov", numObs = nrow(myFADataRaw))
)
m1$S$lbound <- .1
m1 = mxRun(m1)

set.seed(170109)
tmp <- mxTryHard(mxModel(m1, mxCI(c("x1_resid","S[1,1]"), boundAdj=FALSE)), intervals = T)
omxCheckCloseEnough(nrow(tmp$output$confidenceIntervals), 1)

omxCheckCloseEnough(tmp$output$confidenceIntervals[1, c('lbound', 'ubound')],
                    c(.280, .430), .01)

