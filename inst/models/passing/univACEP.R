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


# SCRIPT: univACEP.R
# History:  Sat 26 Sep 2009 16:20:23 BST
#    (tb) Built MZ/DZ models from shared model correctly;  Adapt to data(); Added summary() calls
# TODO: could use omxCheckCloseEnough on the reduced model (tb) 
# OpenMx: http://www.openmx.virginia.com
##########################################

require(OpenMx)
# Prepare Data
data("twinData", package="OpenMx")
selVars <- c('bmi1','bmi2')
aceVars <- c("A1","C1","E1","A2","C2","E2")
mzfData <- subset(twinData, zyg==1, selVars)
dzfData <- subset(twinData, zyg==3, selVars)
cov(mzfData, use="pairwise.complete.obs")
cov(dzfData, use="pairwise.complete.obs")

#Fit ACE Model with RawData and Matrices Input
share <- mxModel("share", type="RAM",
    manifestVars=selVars,
    latentVars=aceVars,
    mxPath(from=aceVars, arrows=2, free=FALSE, values=1),
    mxPath(from="one", to=aceVars, arrows=1, free=FALSE, values=0),
    mxPath(from="one", to=selVars, arrows=1, free=TRUE, values=20, labels= c("mean","mean")),
    mxPath(from="C1",  to="C2",    arrows=2, free=FALSE, values=1),
    mxPath(from=c("A1","C1","E1"), to="bmi1", arrows=1, free=TRUE, values=.6, label=c("a","c","e")),
    mxPath(from=c("A2","C2","E2"), to="bmi2", arrows=1, free=TRUE, values=.6, label=c("a","c","e"))
)

MZ = mxModel(share, name="MZ",
    mxPath(from="A1", to="A2", arrows=2, free=FALSE, values=1),
    mxData(mzfData, type="raw")
)

DZ = mxModel(share, name="DZ",
    mxPath(from="A1", to="A2", arrows=2, free=FALSE, values=.5),
    mxData(dzfData, type="raw") 
)

twinACEModel <- mxModel("twinACE", MZ, DZ,
    mxAlgebra(MZ.objective + DZ.objective, name="twin"), 
    mxFitFunctionAlgebra("twin")
)

#Run ACE model
twinACEFit <- mxRun(twinACEModel)
summary(twinACEFit)

MZc <- twinACEFit$MZ.fitfunction$info$expCov
DZc <- twinACEFit$DZ.fitfunction$info$expCov
M   <- twinACEFit$MZ.fitfunction$info$expMean
A   <- mxEval(a*a, twinACEFit)
C   <- mxEval(c*c, twinACEFit)
E   <- mxEval(e*e, twinACEFit)
V   <- (A+C+E)
a2  <- A/V
c2  <- C/V
e2  <- E/V
ACEest <- rbind(cbind(A,C,E),cbind(a2,c2,e2))
LL_ACE <- mxEval(objective, twinACEFit)

#Mx answers hard-coded
#1: Heterogeneity Model
Mx.A <- 0.6173023
Mx.C <- 5.595822e-14
Mx.E <- 0.1730462
Mx.M <- matrix(c(21.39293, 21.39293),1,2)
Mx.LL_ACE <- 4067.663

#Compare OpenMx results to Mx results (LL: likelihood; EC: expected covariance, EM: expected means)
omxCheckCloseEnough(LL_ACE,Mx.LL_ACE,.001)
omxCheckCloseEnough(A,Mx.A,.001)
omxCheckCloseEnough(C,Mx.C,.001)
omxCheckCloseEnough(E,Mx.E,.001)
omxCheckCloseEnough(M,Mx.M,.001)

#Build and Run AE model
share2 <- mxModel(share,
	mxPath(from=c("C1"), to="bmi1", arrows=1, free=FALSE, values=0, label="c"),
	mxPath(from=c("C2"), to="bmi2", arrows=1, free=FALSE, values=0, label="c")
)

MZ = mxModel(share2, name="MZ",
    mxPath(from="A1", to="A2", arrows=2, free=FALSE, values=1),
    mxData(mzfData, type="raw")
)

DZ = mxModel(share2, name="DZ",
    mxPath(from="A1", to="A2", arrows=2, free=FALSE, values=.5),
    mxData(dzfData, type="raw") 
)

model <- mxModel("twinAE", MZ, DZ,
    mxAlgebra(MZ.objective + DZ.objective, name="twin"), 
    mxFitFunctionAlgebra("twin")
)

fit <- mxRun(model)
summary(fit)

MZc <- fit$MZ.fitfunction$info$expCov
DZc <- fit$DZ.fitfunction$info$expCov
M   <- fit$MZ.fitfunction$info$expMean
A   <- mxEval(a*a, fit)
C   <- mxEval(c*c, fit)
E   <- mxEval(e*e, fit)
V   <- (A + C + E)
a2  <- A / V
c2  <- C / V
e2  <- E / V
AEest <- rbind(cbind(A, C, E),cbind(a2, c2, e2))
LL_AE <- mxEval(objective, fit)

LRT_ACE_AE <- (LL_AE - LL_ACE)

#Print relevant output
ACEest
AEest
LRT_ACE_AE

# NEEDS TEST ANSWERS FROM AE MODEL IN Mx
#Mx answers hard-coded
#1: Heterogeneity Model
# Mx.A <- 0.6173023
# Mx.C <- 5.595822e-14
# Mx.E <- 0.1730462
# Mx.M <- matrix(c(21.39293, 21.39293),1,2)
# Mx.LL_ACE <- 4067.663

#Compare OpenMx results to Mx results (LL: likelihood; EC: expected covariance, EM: expected means)
# omxCheckCloseEnough(LL_ACE,Mx.LL_ACE,.001)
# omxCheckCloseEnough(A,Mx.A,.001)
# omxCheckCloseEnough(C,Mx.C,.001)
# omxCheckCloseEnough(E,Mx.E,.001)
# omxCheckCloseEnough(M,Mx.M,.001)
