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
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2016-02-03
# Filename: univACErSEM.R
# Purpose: Define a behavior genetics single-trait ACE model as a 
#  Relational SEM (rSEM)
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
require(OpenMx)


#------------------------------------------------------------------------------
# Prepare Data

data("twinData", package="OpenMx")
selVars <- c('bmi1','bmi2','zyg')
wideData <- subset(twinData, zyg %in% c(1, 3), selVars)
wideData$rel <- c(1, NA, .5)[wideData$zyg]
wideData$famID <- 1:nrow(wideData)
tallData <- reshape(wideData, varying=c('bmi1', 'bmi2'), v.names='bmi',
		    timevar='twin', times=1:2, idvar='famID', direction='long')
tallData$personID <- 1:nrow(tallData)
tallData$relsqrt <- sqrt(tallData$rel)
tallData$relu <- sqrt(1-tallData$rel)
tallData <- tallData[order(tallData$famID, tallData$twin),
		     c('famID', 'personID', 'twin', 'rel',
		       'relsqrt', 'relu', 'bmi')]
wData <- tallData
bData <- tallData[!duplicated(tallData$famID),
		  c('famID', 'rel', 'relsqrt')]


#------------------------------------------------------------------------------
# Between Model

bModel <- mxModel(
    'between', type="RAM",
    mxData(type="raw", observed=bData, primaryKey="famID"),
    latentVars = c("C", "AC"),
    mxPath("C", arrows=2, values=1, labels="v_C", lbound=1e-6),
    mxPath("AC", arrows=2, values=1, labels="v_A", lbound=1e-6))


#------------------------------------------------------------------------------
# Within Model

wModel <- mxModel(
    'within', type="RAM", bModel,
    mxData(type="raw", observed=wData),
    manifestVars = 'bmi',
    latentVars = c("E", "AU"),
    mxPath(from="one", to="bmi", arrows=1, free=TRUE, values=20, labels="mean"),
    mxPath('E', arrows=2, values=1, labels="v_E", lbound=1e-6),
    mxPath('AU', arrows=2, values=1, labels="v_A", lbound=1e-6),
    mxPath('AU', 'bmi', values=1, labels='data.relu', free=FALSE),
    mxPath('E', 'bmi', free=FALSE, values=1),
    mxPath('between.C', 'bmi', values=1,
	   free=FALSE, joinKey="famID"),
    mxPath('between.AC', 'bmi', values=1, arrows=1, free=FALSE,
	   labels='data.relsqrt', joinKey="famID"))


#------------------------------------------------------------------------------
# Run 'em
wRun <- mxRun(wModel)


#------------------------------------------------------------------------------
# Take a look

summary(wRun)

# Cf. inst/models/passing/univACEP.R

#Mx answers hard-coded
#1: Heterogeneity Model
Mx.A <- 0.6173023
Mx.C <- 5.595822e-14
Mx.E <- 0.1730462
Mx.M <- 21.39293
Mx.LL_ACE <- 4067.663

wparam <- mxEval(rbind(v_A, v_C, v_E, mean), wRun)
mparam <- rbind(Mx.A, Mx.C, Mx.E, Mx.M)
omxCheckCloseEnough(wparam, mparam, .001)

omxCheckCloseEnough(-2*logLik(wRun), Mx.LL_ACE, .001)


#------------------------------------------------------------------------------
# Same model, but with constant between-level transition matrix

bLatent <- c('C', 'AC')
bModel2 <- mxModel(
    'between',
    mxData(type="raw", observed=bData, primaryKey="famID"),
    latentVars = bLatent,
    mxMatrix(name="F", nrow=0, ncol=2, dimnames=list(NULL, bLatent)),
    mxAlgebra(data.rel * v_A, name="rel_v_A"),
    mxMatrix("Symm", name="S", nrow=2, ncol=2, dimnames=list(bLatent,bLatent),
	     free=c(TRUE,FALSE,FALSE), labels=c("v_C", NA, "rel_v_A[1,1]"),
	     values=c(1,0,1), lbound=c(1e-6,NA,1e-6)),
    mxMatrix(name="A", nrow=2, ncol=2, values=0,
	     dimnames=list(bLatent,bLatent)),
    mxFitFunctionML(),
    mxExpectationRAM())

#------------------------------------------------------------------------------
# Within Model

wModel2 <- mxModel(
    'within', type="RAM", bModel2,
    mxData(type="raw", observed=wData),
    manifestVars = 'bmi',
    latentVars = c("E", "AU"),
    mxPath(from="one", to="bmi", arrows=1, free=TRUE,
	   values=20, labels="mean"),
    mxPath('E', arrows=2, values=1, labels="v_E", lbound=1e-6),
    mxPath('AU', arrows=2, values=1, labels="v_A", lbound=1e-6),
    mxPath('AU', 'bmi', values=1, labels='data.relu', free=FALSE),
    mxPath('E', 'bmi', free=FALSE, values=1),
    mxPath('between.C', 'bmi', values=1,
	   free=FALSE, joinKey="famID"),
    mxPath('between.AC', 'bmi', values=1,
	   free=FALSE, joinKey="famID"))

# This isn't a huge speed-up because the per-cluster covariance matrix
# is already small in the version above.
wRun2 <- mxRun(wModel2)

wparam <- mxEval(rbind(v_A, v_C, v_E, mean), wRun2)
mparam <- rbind(Mx.A, Mx.C, Mx.E, Mx.M)
omxCheckCloseEnough(wparam, mparam, .001)

omxCheckCloseEnough(-2*logLik(wRun2), Mx.LL_ACE, .001)

omxCheckCloseEnough(wRun2$expectation$debug$rampartUsage, 867, 1)
