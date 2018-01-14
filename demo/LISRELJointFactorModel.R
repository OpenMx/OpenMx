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


#--------------------------------------
# Fit Joint Ordinal/Continuous Factor Model

require(OpenMx)

# get data (loaded from demo data sets in OpenMx package)
data(jointdata)

# specify ordinal columns as ordered factors
jointdata[,c(2,4,5)] <- mxFactor(jointdata[,c(2,4,5)], 
	levels=list(c(0,1), c(0, 1, 2, 3), c(0, 1, 2)))

loadings <- mxMatrix("Full", 5, 1,
	free=TRUE, values=1, name="L", dimnames=list(names(jointdata), "Factor1"))

resid <- mxMatrix("Diag", 5, 5,
	free=c(TRUE, FALSE, TRUE, FALSE, FALSE), values=.5, name="U")

means <- mxMatrix("Full", 5, 1,
	free=c(TRUE, FALSE, TRUE, FALSE, FALSE), values=0, name="M",
	dimnames=list(names(jointdata), NA))

ident <- mxMatrix("Diag", 1, 1, FALSE, 1, name="I")
zerom <- mxMatrix("Zero", 1, 1, name="Z", dimnames=list("Factor1", NA))

thrFre <- c(TRUE, FALSE, FALSE, rep(TRUE, 5), FALSE)
thrVal <- c(0, NA, NA, -1, 0, 1, -1, 1, NA)
thrLab <- c("z2t1", NA, NA, "z4t1", "z4t2", "z4t3", "z5t1", "z5t2", NA)
thresh <- mxMatrix("Full", 3, 3, free=thrFre, values=thrVal, labels=thrLab,
	name="T", dimnames=list(c(NA, NA, NA), c("z2", "z4", "z5")))


# run factor model
jointModel1 <- mxModel("ContinuousOrdinalData",
	mxData(jointdata[1:90,], "raw"), # limit num rows for example speed
	loadings, resid, means, ident, zerom, thresh,
	mxFitFunctionML(),
	mxExpectationLISREL(LX="L", TX="M", PH="I", KA="Z", TD="U",
		dimnames=names(jointdata),
		thresholds="T",
		threshnames=c("z2", "z4", "z5"))
)

# Run the joint model
jointResults1 <- mxRun(jointModel1, suppressWarnings=TRUE)
