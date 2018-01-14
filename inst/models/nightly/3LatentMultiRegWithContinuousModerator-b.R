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
# Program: 3LatentMultiRegWithModerator100521.R
#  Author: Steven M. Boker
#    Date: Sat May 22 11:13:51 EDT 2010
#
# This program tests variations on a latent variable multiple regression
#    using a standard RAM.
#
# ---------------------------------------------------------------------
# Revision History
#    -- Sat May 22 11:13:51 EDT 2010
#      Created 3LatentMultiRegWithModerator100521.R.
#
# ---------------------------------------------------------------------

# ----------------------------------
# Read libraries and set options.

library(OpenMx)

options(width=100)

# ---------------------------------------------------------------------
# Data for multiple regression of F3 on F1 and F2 with moderator variable Z.

numberSubjects <- 1000
numberIndicators <- 12
numberFactors <- 3

set.seed(10)

fixedBMatrixF <- matrix(c(.4, .2), 2, 1, byrow=TRUE)
randomBMatrixF <- matrix(c(.3, .5), 2, 1, byrow=TRUE)
XMatrixF <- matrix(rnorm(numberSubjects*2, mean=0, sd=1), numberSubjects, 2)
UMatrixF <- matrix(rnorm(numberSubjects*1, mean=0, sd=1), numberSubjects, 1)
Z <- matrix(rnorm(numberSubjects, mean=0, sd=1), nrow=numberSubjects, ncol=2)

XMatrix <- cbind(XMatrixF, XMatrixF %*% fixedBMatrixF + (XMatrixF*Z) %*% randomBMatrixF + UMatrixF)

BMatrix <- matrix(c( 1, .6, .7, .8,  0,  0,  0,  0,  0,  0,  0,  0,
                     0,  0,  0,  0,  1, .5, .6, .7,  0,  0,  0,  0,
                     0,  0,  0,  0,  0,  0,  0,  0,  1, .7, .6, .5), numberFactors, numberIndicators, byrow=TRUE)
UMatrix <- matrix(rnorm(numberSubjects*numberIndicators, mean=0, sd=1), numberSubjects, numberIndicators)
YMatrix <- XMatrix %*% BMatrix + UMatrix

cor(cbind(XMatrix,Z[,1]))

dimnames(YMatrix) <- list(NULL, paste("X", 1:numberIndicators, sep=""))

latentMultiRegModerated1 <- cbind(YMatrix,Z=Z[,1])

round(cor(latentMultiRegModerated1), 3)
round(cov(latentMultiRegModerated1), 3)

latentMultiRegModerated1[,'Z'] <- latentMultiRegModerated1[,'Z'] - mean(latentMultiRegModerated1[,'Z'])

numberFactors <- 3
numberIndicators <- 12
numberModerators <- 1
indicators <- paste("X", 1:numberIndicators, sep="")
moderators <- c("Z")
totalVars <- numberIndicators + numberFactors + numberModerators

# ----------------------------------
# Build an orthogonal simple structure factor model

latents <- paste("F", 1:numberFactors, sep="")

uniqueLabels <- paste("U_", indicators, sep="")
meanLabels <- paste("M_", latents, sep="")
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

threeLatentOrthogonal <- mxModel("threeLatentOrthogonal",
    type="RAM",
    manifestVars=c(indicators),
    latentVars=c(latents,"dummy1"),
    mxPath(from=latents1, to=indicators1, 
           arrows=1, connect="all.pairs",
           free=TRUE, values=.2, 
           labels=loadingLabels1),
    mxPath(from=latents2, to=indicators2, 
           arrows=1, connect="all.pairs",
           free=TRUE, values=.2, 
           labels=loadingLabels2),
    mxPath(from=latents3, to=indicators3, 
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
           free=TRUE, values=.2, 
           labels=uniqueLabels),
    mxPath(from=latents,
           arrows=2, lbound=1e-5,
           free=TRUE, values=.8, 
           labels=factorVarLabels),
    mxPath(from="one", to=indicators, 
           arrows=1, free=FALSE, values=0),
    mxPath(from="one", to=c(latents), 
           arrows=1, free=TRUE, values=.1, 
           labels=meanLabels),
    mxData(observed=latentMultiRegModerated1, type="raw")
    )

# ----------------------------------
# Modify to add in direct paths


threeLatentNoModerator <- mxModel(threeLatentOrthogonal,
    mxPath(from=c("F1","F2"),to="F3",
           arrows=1, lbound=0,
           free=TRUE, values=.2, labels=c("b11", "b12")),
    mxPath(from="F1",to="F2",
           arrows=2, 
           free=TRUE, values=.1, labels=c("cF1F2")),
    name="threeLatentNoModerator"
    )

threeLatentNoModeratorOut <- mxRun(threeLatentNoModerator)
summary(threeLatentNoModeratorOut)
omxCheckCloseEnough(threeLatentNoModeratorOut$output$fit, 37871.72, .1)
