#
#   Copyright 2007-2024 by the individuals mentioned in the source code history
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

# Adapted from script by Timothy C. Bates

library(OpenMx)
# This script doesn't test anything optimization-related, so it doesn't need 
# to be run with anything other than the on-load default:
if(mxOption(NULL, 'Default optimizer') != "SLSQP") stop("SKIP")
mxVersion()
set.seed(241218)
tData = data.frame(matrix(rnorm(300,0,1), nrow=100, ncol=3))
names(tData) = c("t1", "t2", "t3")

svModel = mxModel(
	model="startVals", type="RAM", 
	manifestVars = c("t1", "t2", "t3"),
	mxPath("t1", to = "t2", values = .1, labels = "b"),
	mxPath("t1", to = "t3", values = .2, labels = "b"),
	mxPath(c("t1","t2","t3"), arrows = 2, values = .5, connect = "single"),
	mxPath("one", to = c("t1","t2","t3"), values = 0),
	mxData(tData, type = "raw")
)
omxCheckError(
	mxRun(svModel),
	"The free parameter 'b' has been assigned multiple starting values! See matrix 'A' at location (3, 1) and matrix 'A' at location (2, 1) If you want to randomly select one of these values, call model <- omxAssignFirstParameters(model) before running again."
)
svModel <- omxAssignFirstParameters(svModel)
mxRun(svModel) #<--Should run without error.

lbModel = mxModel(
	model="Lbounds", type="RAM", 
	manifestVars = c("t1", "t2", "t3"),
	mxPath("t1", to = "t2", values = .1, labels = "b", lbound = .2),
	mxPath("t1", to = "t3", values = .1, labels = "b", lbound = .1),
	mxPath(c("t1","t2","t3"), arrows = 2, values = .5, connect = "single"),
	mxPath("one", to = c("t1","t2","t3"), values = 0),
	mxData(tData, type = "raw")
)
omxCheckError(
	mxRun(lbModel),
	"The free parameter 'b' has been assigned multiple lower bounds! See matrix 'A' at location (3, 1) and matrix 'A' at location (2, 1) If you want to randomly select one of these lbounds, call model <- omxAssignFirstParameters(model) before running again."
)
lbModel <- omxAssignFirstParameters(lbModel)
mxRun(lbModel) #<--Should run without error.

ubModel = mxModel(
	model="Ubounds", type="RAM", 
	manifestVars = c("t1", "t2", "t3"),
	mxPath("t1", to = "t2", values = .1, labels = "b", ubound = .2),
	mxPath("t1", to = "t3", values = .1, labels = "b", ubound = .1),
	mxPath(c("t1","t2","t3"), arrows = 2, values = .5, connect = "single"),
	mxPath("one", to = c("t1","t2","t3"), values = 0),
	mxData(tData, type = "raw")
)
omxCheckError(
	mxRun(ubModel),
	"The free parameter 'b' has been assigned multiple upper bounds! See matrix 'A' at location (3, 1) and matrix 'A' at location (2, 1) If you want to randomly select one of these ubounds, call model <- omxAssignFirstParameters(model) before running again."
)
ubModel <- omxAssignFirstParameters(ubModel)
mxRun(ubModel) #<--Should run without error.

# Test with more than two offending instances: ####

svModel = mxModel(
	model="startVals", type="RAM", 
	manifestVars = c("t1", "t2", "t3"),
	mxPath("t1", to = "t2", values = .1),
	mxPath("t1", to = "t3", values = .2),
	#mxPath("t1", to = "t3", values = .3, labels = "b"),
	mxPath(c("t1","t2","t3"), arrows = 2, values = .5, connect = "single"),
	mxPath("one", to = c("t1","t2","t3"), values = c(0,0.1,0.2), labels="b"),
	mxData(tData, type = "raw")
)
omxCheckError(
	mxRun(svModel),
	"The free parameter 'b' has been assigned multiple starting values! See matrix 'M' at location (1, 2) and matrix 'M' at location (1, 1) If you want to randomly select one of these values, call model <- omxAssignFirstParameters(model) before running again."
)
svModel <- omxAssignFirstParameters(svModel)
mxRun(svModel) #<--Should run without error.

lbModel = mxModel(
	model="Lbounds", type="RAM", 
	manifestVars = c("t1", "t2", "t3"),
	mxPath("t1", to = "t2", values = .1),
	mxPath("t1", to = "t3", values = .1),
	#mxPath("t1", to = "t3", values = .1, labels = "b", lbound = .3),
	mxPath(c("t1","t2","t3"), arrows = 2, values = .5, connect = "single"),
	mxPath("one", to = c("t1","t2","t3"), values = 0, lbound=c(0,0.1,0.2), labels="b"),
	mxData(tData, type = "raw")
)
omxCheckError(
	mxRun(lbModel),
	"The free parameter 'b' has been assigned multiple lower bounds! See matrix 'M' at location (1, 2) and matrix 'M' at location (1, 1) If you want to randomly select one of these lbounds, call model <- omxAssignFirstParameters(model) before running again."
)
lbModel <- omxAssignFirstParameters(lbModel)
mxRun(lbModel) #<--Should run without error.

ubModel = mxModel(
	model="Ubounds", type="RAM", 
	manifestVars = c("t1", "t2", "t3"),
	mxPath("t1", to = "t2", values = .1, labels = "b", ubound = .2),
	mxPath("t1", to = "t3", values = .1, labels = "b", ubound = .1),
	#mxPath("t2", to = "t1", values = .1, labels = "b", ubound = .3),
	mxPath(c("t1","t2","t3"), arrows = 2, values = .5, connect = "single"),
	mxPath("one", to = c("t1","t2","t3"), values = 0, lbound=c(1.0,1.1,1.2), labels="b"),
	mxData(tData, type = "raw")
)
omxCheckError(
	mxRun(ubModel),
	"The free parameter 'b' has been assigned multiple upper bounds! See matrix 'A' at location (3, 1) and matrix 'A' at location (2, 1) If you want to randomly select one of these ubounds, call model <- omxAssignFirstParameters(model) before running again."
)
ubModel <- omxAssignFirstParameters(ubModel)
mxRun(ubModel)

# Test multiple things wrong: ####

tModel = mxModel(
	model="startVals", type="RAM", 
	manifestVars = c("t1", "t2", "t3"),
	mxPath("t1", to = "t2", values = .1, labels = "b", lbound=0, ubound=1),
	mxPath("t1", to = "t3", values = .2, labels = "b", lbound=0.1, ubound=1.1),
	mxPath(c("t1","t2","t3"), arrows = 2, values = .5, connect = "single"),
	mxPath("one", to = c("t1","t2","t3"), values = 0),
	mxData(tData, type = "raw")
)
omxCheckError(
	mxRun(tModel),
	"The free parameter 'b' has been assigned multiple starting values! See matrix 'A' at location (3, 1) and matrix 'A' at location (2, 1) If you want to randomly select one of these values, call model <- omxAssignFirstParameters(model) before running again."
)
tModel <- omxAssignFirstParameters(tModel)
mxRun(tModel) #<--Should run without error.
