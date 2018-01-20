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
# Program: RAMInverseTest.R
#  Author: Michael D. Hunter
#    Date: 2012.03.23
#
# This program tests the RAM inverse speedup using a latent variable
#  regression example.
#
# ---------------------------------------------------------------------
# Revision History
#  Fri Mar 23 22:45:52 Central Daylight Time 2012 -- Michael Hunter Created RAMInverseTest.R
#   from Steven Boker's test IntroSEM-ThreeLatentMultipleRegTest1.R
#  Sat Mar 24 01:40:52 Central Daylight Time 2012 -- Michael Hunter changed model to use ML instead of FIML
#  
#
# ---------------------------------------------------------------------


#------------------------------------------------------------------------
#------------------------------------------------------------------------
#------------------------------------------------------------------------
# BEGIN COPIED SECTION FROM model/passing/IntroSEM-ThreeLatentMultipleRegTest1.R


# ----------------------------------
# Read libraries and set options.

require(OpenMx)


# ----------------------------------
# Read the data and print descriptive statistics.

data(latentMultipleRegExample1)

numberFactors <- 3
indicators <- names(latentMultipleRegExample1)
numberIndicators <- length(indicators)
totalVars <- numberIndicators + numberFactors


# ----------------------------------
# Build an Old-style RAM OpenMx single factor FIML model with fixed variance

latents <- paste("F", 1:numberFactors, sep="")
loadingLabels <- paste("b_F", rep(1:numberFactors, each=numberIndicators), 
                              rep(indicators, numberFactors), sep="") 
loadingLabels

uniqueLabels <- paste("U_", indicators, sep="")
meanLabels <- paste("M_", indicators, sep="")
factorVarLabels <- paste("Var_", latents, sep="")

latents1 <- latents[1]
indicators1 <- indicators[1:4]
loadingLabels1 <- paste("b_F1", indicators[1:4], sep="") 
latents2 <- latents[2]
indicators2 <- indicators[5:8]
loadingLabels2 <- paste("b_F2", indicators[5:8], sep="") 
latents3 <- latents[3]
indicators3 <- indicators[9:12]
loadingLabels3 <- paste("b_F3", indicators[9:12], sep="") 

threeLatentOrthoRaw1 <- mxModel("threeLatentOrthogonal",
    type="RAM",
    manifestVars=indicators,
    latentVars=latents,
    mxPath(from=latents1, to=indicators1, 
#           arrows=1, all=TRUE,
           arrows=1, connect="all.pairs",
           free=TRUE, values=.2, 
           labels=loadingLabels1),
    mxPath(from=latents2, to=indicators2, 
#           arrows=1, all=TRUE, 
           arrows=1, connect="all.pairs", 
           free=TRUE, values=.2, 
           labels=loadingLabels2),
    mxPath(from=latents3, to=indicators3, 
#           arrows=1, all=TRUE,
           arrows=1, connect="all.pairs", 
           free=TRUE, values=.2, 
           labels=loadingLabels3),
    mxPath(from=latents1, to=indicators1[1], 
           arrows=1, 
           free=FALSE, values=1),
    mxPath(from=latents2, to=indicators2[1], 
           arrows=1, 
           free=FALSE, values=1),
    mxPath(from=latents3, to=indicators3[1], 
           arrows=1, 
           free=FALSE, values=1),
    mxPath(from=indicators, 
           arrows=2, 
           free=TRUE, values=.8, 
           labels=uniqueLabels),
    mxPath(from=latents,
           arrows=2, 
           free=TRUE, values=.8, 
           labels=factorVarLabels),
    mxPath(from="one", to=indicators, 
           arrows=1, free=TRUE, values=.1, 
           labels=meanLabels),
    mxData(observed=cov(latentMultipleRegExample1), type="cov", numObs=nrow(latentMultipleRegExample1), means=colMeans(latentMultipleRegExample1))
    )

threeLatentOrthoRaw1Out <- mxRun(threeLatentOrthoRaw1, suppressWarnings=TRUE)

summary(threeLatentOrthoRaw1Out)

# model 2
threeLatentObliqueRaw1 <- mxModel(threeLatentOrthoRaw1,
#    mxPath(from=latents,to=latents,all=TRUE,
    mxPath(from=latents,to=latents,connect="unique.pairs",
           arrows=2, 
           free=TRUE, values=.3),
    mxPath(from=latents,
           arrows=2, 
           free=TRUE, values=.8, 
           labels=factorVarLabels),
    name="threeLatentOblique"
    )

threeLatentObliqueRaw1Out <- mxRun(threeLatentObliqueRaw1, suppressWarnings=TRUE)

summary(threeLatentObliqueRaw1Out)

# model 3
threeLatentMultipleReg1 <- mxModel(threeLatentOrthoRaw1,
    mxPath(from="F1",to="F2",
           arrows=2, 
           free=TRUE, values=.3),
    mxPath(from=c("F1","F2"), to="F3",
           arrows=1, 
           free=TRUE, values=.2, 
           labels=c("b1", "b2")),
    name="threeLatentMultipleReg"
    )

threeLatentMultipleReg1Out <- mxRun(threeLatentMultipleReg1, suppressWarnings=TRUE)

# END COPIED SECTION FROM models/passing/IntroSEM-ThreeLatentMultipleRegTest1.R
#------------------------------------------------------------------------
#------------------------------------------------------------------------
#------------------------------------------------------------------------


#------------------------------------------------------------------------
# Note that the A matrix is nilpotent of nontrivial order, 3

AM <- threeLatentMultipleReg1$A$values
AM
AM %*% AM
AM %*% AM %*% AM


#------------------------------------------------------------------------
# Turn off the RAM inverse speed up

threeLatentMultipleReg1B <- mxOption(threeLatentMultipleReg1, "RAM Inverse Optimization", "No")
brun <- mxRun(threeLatentMultipleReg1B) # This causes an error when the fast inverse is different from the slow


#------------------------------------------------------------------------
# Compare the run times of the same model using the RAM speed up or not
summary(threeLatentMultipleReg1Out)$wallTime
summary(brun)$wallTime


#------------------------------------------------------------------------
# Compare the estimated parameters of the same model using the speed up or not

omxCheckCloseEnough(threeLatentMultipleReg1Out$output$estimate, brun$output$estimate, epsilon=0.001)


#------------------------------------------------------------------------
# End
