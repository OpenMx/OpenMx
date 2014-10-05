#
#   Copyright 2007-2012 The OpenMx Project
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

# -----------------------------------------------------------------------------

require(OpenMx)
	
data(myFADataRaw)

oneFactorOrd <- myFADataRaw[,c("z1", "z2", "z3")]

oneFactorOrd$z1 <- mxFactor(oneFactorOrd$z1, levels=c(0, 1))
oneFactorOrd$z2 <- mxFactor(oneFactorOrd$z2, levels=c(0, 1))
oneFactorOrd$z3 <- mxFactor(oneFactorOrd$z3, levels=c(0, 1, 2))

oneFactorModel <- mxModel("Common Factor Model Path Specification", type = "RAM",
	mxData(observed=oneFactorOrd, type = "raw"),
	manifestVars=c("z1","z2","z3"),
	latentVars="F1",
	# residual variances
	mxPath(from=c("z1","z2","z3"), arrows=2, free=F, values=c(1,1,1), labels=c("e1","e2","e3")), # latent variance
	mxPath(from="F1", arrows=2, free=T, values=1, labels ="varF1"), # factor loadings
	mxPath(from="F1", to=c("z1","z2","z3"), arrows=1, free=c(F,T,T), values=c(1,1,1), labels=c("l1","l2","l3")), # means
	mxPath(from="one", to=c("z1","z2","z3","F1"), arrows=1, free=F, values=0, labels=paste("meanz", c(1:3, "F"))), # thresholds
	mxThreshold(vars=c("z1", "z2", "z3"), nThresh=c(1,1,2), free=T, values=c(-1, 0, -.5, 1.2))
) # close model

oneFactorResults <- mxRun(oneFactorModel)

# joint ordinal-continuous model

oneFactorJoint <- myFADataRaw[,c("x1", "x2", "x3", "z1", "z2", "z3")]
	
oneFactorJoint$z1 <- mxFactor(oneFactorOrd$z1, levels=c(0, 1))
oneFactorJoint$z2 <- mxFactor(oneFactorOrd$z2, levels=c(0, 1))
oneFactorJoint$z3 <- mxFactor(oneFactorOrd$z3, levels=c(0, 1, 2))

oneFactorJointModel <- mxModel("Common Factor Model Path Specification", type="RAM",
	mxData(observed=oneFactorJoint, type="raw"),
	manifestVars=c("x1", "x2", "x3", "z1","z2","z3"),
	latentVars="F1",
	# residual variances
	mxPath(from=c("x1", "x2", "x3", "z1","z2","z3"), arrows=2, 
		   free=c(T, T, T, F, F, F), values=1, 
		   labels=c("e1","e2","e3","e4","e5","e6")), # latent variance
	mxPath(from="F1", arrows=2, free=F, values=1, labels ="varF1"), # factor loadings
	mxPath(from="F1", to=c("x1", "x2", "x3", "z1","z2","z3"), arrows=1, free=T, values=1, 
	       labels=c("l1","l2","l3","l4","l5","l6")), # means
	mxPath(from="one", to=c("x1", "x2", "x3","z1","z2","z3","F1"), arrows=1, 
	       free=c(T,T,T,F,F,F,F), values=0, 
		   labels=c("meanx1","meanx2","meanx3","meanz1","meanz2","meanz3","meanF")), # thresholds
	mxThreshold(vars=c("z1", "z2", "z3"), nThresh=c(1,1,2), free=T, values=c(-1, 0, -.5, 1.2))
) # close model

oneFactorJointModelFit <- mxRun(oneFactorJointModel)
