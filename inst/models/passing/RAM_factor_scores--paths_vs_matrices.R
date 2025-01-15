#
#   Copyright 2007-2025 by the individuals mentioned in the source code history
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

# This test script verifies that `mxFactorScores()` gives substantially identical results
# for a RAM model specified using paths vis-a-vis matrices.

require(OpenMx)
data(demoOneFactor)
manifests <- names(demoOneFactor)
latents <- c("G")
fm1 <- mxModel(
	"OneFactor",
	type="RAM",
	manifestVars = manifests,
	latentVars = latents,
	mxPath(from=latents, to=manifests,values=0.8),
	mxPath(from=manifests, arrows=2,values=1),
	mxPath(from=latents, arrows=2,
				 free=FALSE, values=1.0),
	mxPath(from="one",to=manifests,values=0.1),
	mxData(demoOneFactor, type="raw")
)
fmf1 <- mxRun(fm1)
fsml1 <- mxFactorScores(fmf1,"ML")
fswml1 <- mxFactorScores(fmf1,"WeightedML")
fsr1 <- mxFactorScores(fmf1,"Regression")

fm2 <- mxModel(
	"OneFactor",
	mxMatrix(type="Full",nrow=6,ncol=6,free=fm1$A$free,values=fm1$A$values,name="A"),
	mxMatrix(type="Symm",nrow=6,ncol=6,free=fm1$S$free,values=fm1$S$values,name="S"),
	mxMatrix(type="Full",nrow=5,ncol=6,free=F,values=fm1$F$values,name="F"),
	mxMatrix(type="Full",nrow=1,ncol=6,free=fm1$M$free,values=fm1$M$values,name="M"),
	mxExpectationRAM(A="A",S="S",F="F",M="M",dimnames=c(manifests,latents)),
	mxFitFunctionML(),
	mxData(demoOneFactor, type="raw")
)
fmf2 <- mxRun(fm2)
fsml2 <- mxFactorScores(fmf2,"ML")
fswml2 <- mxFactorScores(fmf2,"WeightedML")
fsr2 <- mxFactorScores(fmf2,"Regression")

omxCheckCloseEnough(cor(fsml1[,1,1],fsml2[,1,1]),1,1e-8)
omxCheckCloseEnough(cor(fsml1[,1,2],fsml2[,1,2]),1,1e-8)
omxCheckCloseEnough(cor(fswml1[,1,1],fswml2[,1,1]),1,1e-8)
omxCheckCloseEnough(cor(fswml1[,1,2],fswml2[,1,2]),1,1e-8)
omxCheckCloseEnough(cor(fsr1[,1,1],fsr2[,1,1]),1,1e-8)
omxCheckTrue(all(summary(fsr1[,1,2]-fsr2[,1,2])==0))
