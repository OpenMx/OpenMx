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
data(demoOneFactor)

# for on.exit()
run <- function() {
	manifests <- names(demoOneFactor)
	latents <- c("G1")
	model1 <- mxModel("One Factor", type="RAM",
		manifestVars = manifests,
		latentVars = latents,
		mxPath(from=latents, to=manifests),
		mxPath(from=manifests, arrows=2),
		mxPath(from=latents, arrows=2, free=F, values=1.0),
		mxData(cov(demoOneFactor), type="cov",numObs=500))
	generated <- omxGraphviz(model1, dotFilename=NULL)
	referenceFile <- file("data/one-factor-reference.dot", open="r")
	on.exit(c(close(referenceFile)))
	reference <- readLines(referenceFile)
	omxCheckTrue(identical(generated, paste(paste(reference, sep="", collapse="\n"), sep="", "\n")))
}

run()
