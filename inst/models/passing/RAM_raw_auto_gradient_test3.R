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
	mxPath(from="one",to=manifests,values=0.1),
	mxData(demoOneFactor, type="raw")
)

# First, make sure results match, with versus without analytic derivatives, with complete data ####

mxOption(NULL,"Analytic gradients","Yes")
fmf1 <- mxRun(factorModel)
summary(fmf1)
mxOption(NULL,"Analytic gradients","No")
fmf2 <- mxRun(factorModel)
summary(fmf2)
omxCheckCloseEnough(coef(fmf1),coef(fmf2),5e-6)
omxCheckCloseEnough(fmf1$output$standardErrors,fmf2$output$standardErrors,5e-6)
omxCheckCloseEnough(fmf1$output$fit,fmf2$output$fit,1e-8)

# Now, make sure results match, with versus without analytic derivatives, with missing data ####

demoOneFactor[1:10,1] <- NA
demoOneFactor[11:20,2] <- NA
demoOneFactor[21:30,3] <- NA
demoOneFactor[31:40,4] <- NA
demoOneFactor[41:50,5] <- NA
factorModel <- mxModel(
	"One Factor",
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
mxOption(NULL,"Analytic gradients","Yes")
fmf3 <- mxRun(factorModel)
summary(fmf3)
mxOption(NULL,"Analytic gradients","No")
fmf4 <- mxRun(factorModel)
summary(fmf4)
omxCheckCloseEnough(coef(fmf3),coef(fmf4),5e-6)
omxCheckCloseEnough(fmf3$output$standardErrors,fmf4$output$standardErrors,5e-6)
omxCheckCloseEnough(fmf3$output$fit,fmf4$output$fit,1e-8)

# Now, make sure results match, with versus without analytic derivatives, with missing-data patterns too small to ####
# constitute new sufficientSets.                                                                                  ####

data(demoOneFactor)
demoOneFactor[1,1] <- NA
demoOneFactor[2,2] <- NA
demoOneFactor[3,3] <- NA
demoOneFactor[4,4] <- NA
demoOneFactor[41:50,5] <- NA
#plan <- mxComputeSequence(list(mxComputeOnce("fitfunction",c("fit","gradient")),mxComputeReportDeriv(),mxComputeReportExpectation()))
factorModel <- mxModel(
	"One Factor",
	type="RAM",
#	plan,
	manifestVars = manifests,
	latentVars = latents,
	mxPath(from=latents, to=manifests,values=0.8),
	mxPath(from=manifests, arrows=2,values=1),
	mxPath(from=latents, arrows=2,
				 free=FALSE, values=1.0),
	mxPath(from="one",to=manifests,values=0.1),
	mxData(demoOneFactor, type="raw")
)
mxOption(NULL,"Analytic gradients","Yes")
fmf5 <- mxRun(factorModel)
summary(fmf5)
mxOption(NULL,"Analytic gradients","No")
fmf6 <- mxRun(factorModel)
summary(fmf6)
omxCheckCloseEnough(coef(fmf5),coef(fmf6),5e-6)
omxCheckCloseEnough(fmf5$output$standardErrors,fmf6$output$standardErrors,5e-6)
omxCheckCloseEnough(fmf5$output$fit,fmf6$output$fit,1e-8)
