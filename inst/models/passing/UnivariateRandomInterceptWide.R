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

# ---------------------------------------------------------------------
# Program: UniRandomIntTest-120815.R
#  Author: Steve Boker
#    Date: Wed Aug 15 10:50:12 CEST 2012
#
# This program simulates some univariate multilevel data with random
# intercepts only, fits it with lme(), fits a naive wide format
# multilevel OpenMx model and checks the results
#
# ---------------------------------------------------------------------
# Revision History
#   Steve Boker -- Wed Aug 15 10:50:14 CEST 2012
#      Created UniRandomIntTest-120815.R
#
# ---------------------------------------------------------------------

# ----------------------------------
# Read libraries and set options.

options(width=110)
library(nlme)
library(OpenMx)

# ----------------------------------
# Set constants.

sdLevelOneE <- sqrt(.2)
sdIntercepts <- sqrt(.5)
sdX <- sqrt(1)

N <- 400    # number of participants
P <- 100    # number of observations per participant
b0 <- .5    # Fixed effect intercept
b1 <- .8    # Fixed effect slope

set.seed(1)

# ----------------------------------
# Simulate the data.

X <- rnorm(N*P, 0, sd=sdX)
ID <- rep(1:N, each=P)
b0i <- b0 + rnorm(N, 0, sd=sdIntercepts)
Y <- rep(b0i, each=P) + b1*X + rnorm(N*P, 0, sd=sdLevelOneE)

SimUniRandomIntFrame <- data.frame(ID, X, Y)

# ----------------------------------
# Test with lme().

lmeOut <- summary(lme(Y ~ X, random= list(~ 1 | ID),
		      data=SimUniRandomIntFrame))

# For lme4, use:
# lmerOut <- lmer(Y ~ X + (1 | ID), data=SimUniRandomIntFrame)

# ----------------------------------
# Set constants.

theIDs <- unique(SimUniRandomIntFrame$ID)
totalN <- length(theIDs)
totalVars <- 2

maxP <- 0
for (tID in theIDs) {
    tmask <- SimUniRandomIntFrame$ID==tID
    tLen <- length(SimUniRandomIntFrame$ID[tmask])
    if (tLen > maxP) 
        maxP <- tLen
}

# ----------------------------------
# Wide-format the data frame from tall format.

wideMatrix <- matrix(NA, nrow=totalN, ncol=1 + (maxP*totalVars))
colnames(wideMatrix) <- c("ID", paste("Y",1:maxP, sep=""),
			  paste("X",1:maxP, sep=""))
i <- 1
for (tID in theIDs) {
    wideMatrix[i, 1] <- tID
    tY <- SimUniRandomIntFrame$Y[SimUniRandomIntFrame$ID==tID]
    wideMatrix[i, 2:(length(tY)+1)] <- tY
    tX <- SimUniRandomIntFrame$X[SimUniRandomIntFrame$ID==tID]
    wideMatrix[i, (2+maxP):(length(tY)+1+maxP)] <- tX
    i <- i + 1
}
wideFrame <- data.frame(wideMatrix)

manifestNames <- colnames(wideFrame)[2:dim(wideFrame)[2]]
xNames <- paste("X",1:maxP, sep="")
yNames <- paste("Y",1:maxP, sep="")
latentNames <- c("b0i")

# ----------------------------------
# Build the OpenMx wide model.

OpenMxModelUniRandomIntModel1 <-
  mxModel("OpenMxModelUniRandomIntModel1",
	type="RAM", 
	manifestVars=manifestNames,
    latentVars=latentNames,
    mxPath(from=xNames, to=yNames, connect="single", arrows=1,
	   free=TRUE, values=.2, labels="b1"),
    mxPath(from=xNames, to=xNames, connect="single", arrows=2,
	   free=TRUE, values=.8, labels="vX"),
    mxPath(from=yNames, to=yNames, connect="single", arrows=2,
	   free=TRUE, values=.8, labels="eY"),
    mxPath(from=latentNames, to=yNames, arrows=1, free=FALSE, values=1),
    mxPath(from=latentNames, to=latentNames, connect="single", arrows=2,
	   free=TRUE, values=.8, labels="vb0i"),
    mxPath(from="one", to=c(xNames), arrows=1,
	   free=TRUE, values=1, labels="mX"),
    mxPath(from="one", to=c(latentNames), arrows=1,
	   free=TRUE, values=1, labels="mb0i"),
    mxData(observed=wideFrame, type="raw")
)

# ----------------------------------
# Fit the model and examine the summary results.

omxFit <- mxRun(OpenMxModelUniRandomIntModel1)

summary(omxFit)

omxCheckCloseEnough(lmeOut$coefficients$fixed[1],
		    mxEval(mb0i, model=omxFit), 0.001)

omxCheckCloseEnough(lmeOut$coefficients$fixed[2],
		    mxEval(b1, model=omxFit), 0.001)

omxCheckCloseEnough(lmeOut$sigma,
		    mxEval(sqrt(eY), model=omxFit), 0.001)

omxCheckCloseEnough(sd(c(lmeOut$coefficients$random$ID)),
		    mxEval(sqrt(vb0i), model=omxFit), 0.001)

if (0) {
	omxCheckCloseEnough(lmeOut$coefficients$fixed,
			    fixef(lmerOut), 1e-4)
  omxCheckCloseEnough(lmeOut$sigma, sigma(lmerOut), 1e-4)
  omxCheckCloseEnough(c(lmeOut$coefficients$random$ID),
		      ranef(lmerOut)$ID[[1]], 1e-4)
}
