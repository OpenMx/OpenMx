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

#--------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2013.09.24
# Filename: StateSpaceMissingData.R
# Purpose: Test the state space expectation with missing data.
#  I compare an SSM with no missing data to one with a little mcar missing data.
#            an SSM as a factor model to a LISREL factor model with no missing data.
#            an SSM as a factor model to a LISREL factor model with mcar missing data.
#--------------------------------------------------------------------


#--------------------------------------------------------------------
# Revision History
# Tue Sep 24 12:37:20 Central Daylight Time 2013 -- Michael Hunter created file
# Mon 17 Feb 2014 18:57:17 Central Standard Time -- Michael Hunter added dimnames and variable re-ordering



#------------------------------------------------------------------------------
require(OpenMx)
data(demoOneFactor)
demoOneFactorMiss <- as.matrix(demoOneFactor)
nvar <- ncol(demoOneFactor)
varnames <- colnames(demoOneFactor)
nmiss <- 20

set.seed(23)

missmat <- matrix(c(
	sample(1:nrow(demoOneFactor), size=nmiss),
	sample(1:nvar, size=nmiss, replace=TRUE)), ncol=2, nrow=nmiss)
demoOneFactorMiss[missmat] <- NA

# Try case of plenty of missing data in last row
demoOneFactorMiss[cbind(rep(nrow(demoOneFactorMiss), 3), c(1,3, 4))] <- NA


#------------------------------------------------------------------------------
# Compare an SSM with no missing data to one with a little mcar missing data.

ssModel <- mxModel(name="State Space Manual Example",
	mxMatrix("Full", 1, 1, TRUE, .3, name="A"),
	mxMatrix("Zero", 1, 1, name="B"),
	mxMatrix("Full", nvar, 1, TRUE, .6, name="C", dimnames=list(varnames, "F1")), #Note: dimnames map rows of C matrix to columns of data!
	mxMatrix("Zero", nvar, 1, name="D"),
	mxMatrix("Diag", 1, 1, FALSE, 1, name="Q"),
	mxMatrix("Diag", nvar, nvar, TRUE, .2, name="R"),
	mxMatrix("Zero", 1, 1, name="x0"),
	mxMatrix("Diag", 1, 1, FALSE, 1, name="P0"),
	mxMatrix("Zero", 1, 1, name="u"),
	mxData(observed=demoOneFactor[,sample(size=5, x=names(demoOneFactor), replace=FALSE)], type="raw"),
	mxExpectationStateSpace("A", "B", "C", "D", "Q", "R", "x0", "P0", "u"),
	mxFitFunctionML()
)
ssRun <- mxRun(ssModel)


ssMiss <- mxModel(ssModel, name="With Missing",
	mxData(observed=demoOneFactorMiss, type="raw")
	)

ssMissRun <- mxRun(ssMiss)

ssNoMissParam <- summary(ssRun)$parameters[, c(5, 6)]
ssMissParam <- summary(ssMissRun)$parameters[, c(5, 6)]

# Check if missing data estimates are close to complete data ones
for(i in 1:nrow(ssNoMissParam)){
omxCheckCloseEnough(ssNoMissParam[i,], ssMissParam[i,], epsilon=0.012)
}

# Note: Even though the data columns are in different orders
#  the parameter estimates are the same because the dimnames of C
#  effectively re-arrange the columns of the data.
head(ssRun$data$observed)
head(ssMissRun$data$observed)



#------------------------------------------------------------------------------
# Compare an SSM as a factor model to a LISREL factor model with no missing data.

ssFactor <- mxModel(ssModel, name="As Factor",
	mxMatrix("Full", 1, 1, FALSE, 0, name="A")
	)

ssFactorRun <- mxRun(ssFactor)

liFactor <- mxModel(ssFactor, name="LISREL Factor",
	mxMatrix("Full", 1, 1, FALSE, 0, name="KA", dimnames=list("F1", "F1")),
	mxMatrix("Full", nvar, 1, TRUE, 0, name="TX", dimnames=list(varnames, NA)),
	mxExpectationLISREL(LX="C", PH="Q", TD="R", TX="TX", KA="KA")
	)
liRun <- mxRun(liFactor)
# summary(liRun)


ssFactorNoMissParam <- summary(ssFactorRun)$parameters[1:10, c(5, 6)]
liFactorNoMissParam <- summary(liRun)$parameters[1:10, c(5, 6)]

omxCheckCloseEnough(ssFactorNoMissParam, liFactorNoMissParam, epsilon=0.01)


#------------------------------------------------------------------------------
# Compare an SSM as a factor model to a LISREL factor model with mcar missing data.

ssFactorMiss <- mxModel(ssFactor, name="As Factor with Missing",
	mxData(observed=demoOneFactorMiss, type="raw")
	)

liFactorMiss <- mxModel(liFactor, name="LISREL Factor with Missing",
	mxData(observed=demoOneFactorMiss, type="raw")
	)

ssFactorMissRun <- mxRun(ssFactorMiss)
liMissRun <- mxRun(liFactorMiss)


ssFactorMissParam <- summary(ssFactorMissRun)$parameters[1:10, c(5, 6)]
liFactorMissParam <- summary(liMissRun)$parameters[1:10, c(5, 6)]

omxCheckCloseEnough(ssFactorMissParam, liFactorMissParam, epsilon=0.01)


#------------------------------------------------------------------------------
# Ensure we ignore all-missing rows

demoOneFactorMiss[10,] <- NA

ssFactorAllMiss <- mxModel(
	ssModel, name="As Factor",
	mxMatrix("Full", 1, 1, FALSE, 0, name="A"),
	mxData(observed=demoOneFactorMiss, type="raw"))

ssFactorAllMissRun <- mxRun(ssFactorAllMiss)

param <- summary(ssFactorAllMissRun)$parameters[1:10, c(5, 6)]
omxCheckCloseEnough(param, liFactorMissParam, epsilon=0.01)
