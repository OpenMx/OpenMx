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

require(OpenMx)
require(nlme)

# Multilevel Long Format Test
# Author: Steve Boker
# Date: Sun Nov 29 14:06:07 EST 2009


# This script is used to test the multilevel long format functionality using definition variables
#    as indices.
totalOccasions <- 100
totalSubjects <- 10
set.seed(42) # repeatibility
tID <- rep(1:totalSubjects, each=totalOccasions)
trueX <- rep(rnorm(totalOccasions, mean=0, sd=2), each=totalSubjects) + rnorm(totalOccasions*totalSubjects, mean=0, sd=.2)
trueB <- rep(rnorm(totalSubjects, mean=.8, sd=.3), each=totalOccasions)
tDataFrame <- data.frame(ID=tID, X=trueX, Y=trueB*trueX + rnorm(totalOccasions*totalSubjects,mean=0, sd=.1),trueB=trueB)
summary(tDataFrame)

manifestVars <- c("X", "Y")
numSubjects <- length(unique(tDataFrame$ID))


multilevelModel2 <- mxModel("Multilevel_2",
    mxMatrix("Full", nrow=numSubjects, ncol=2,
        values=c(.2,0),
        free=c(TRUE, TRUE), 
        name="Rand",
        byrow=TRUE
    ),
    mxMatrix("Full", 2, 2, 
        labels=c(NA,  NA,
                "randrow[1,1]", NA), 
        free=FALSE, 
        name="A", 
        byrow=TRUE
    ),
    mxMatrix("Symm", 2, 2,
        values=c(.9,0,.9),
        free=c(T,
               F, T), 
        labels=c("varX",
                 NA, "varY"), 
        name="S", 
        byrow=TRUE
    ),
	mxMatrix("Full", 2, 2, 
	         values=c(1,0,
		              0,1),
		     free=FALSE, 
		     byrow=TRUE, name="F"),
    mxMatrix("Iden", 2, name="I"),
    mxAlgebra(F %*% solve(I-A) %*% S %*% t(solve(I-A)) %*% t(F), 
        name="R", 
        dimnames = list(manifestVars, manifestVars)
    ),
    mxMatrix("Full", nrow=1, ncol=length(manifestVars),
        values=0,
        free=FALSE,
        labels=c(NA,"randrow[1,2]"),
        dimnames=list(NULL, manifestVars),
        name="M"
    ),
	mxAlgebra(Rand[data.ID,], name="randrow"),
    mxFIMLObjective(covariance="R", means="M"),
    mxData(tDataFrame, type="raw")
)

# ----------------------------------
# Fit the  model and examine the summary results.

multilevelModel2Fit <- mxRun(multilevelModel2)

summary(multilevelModel2Fit)

lmeOut <- lme(Y~X, random= ~ X | ID, data=tDataFrame)

cbind(multilevelModel2Fit@output$estimate[1:numSubjects], 
      lmeOut$coef$random$ID[,2] + lmeOut$coef$fixed[2],
      trueB[seq(1,totalOccasions*(totalSubjects), by=totalOccasions)])

mean(multilevelModel2Fit@output$estimate[1:numSubjects])

omxCheckCloseEnough(mean(multilevelModel2Fit@output$estimate[1:numSubjects]), 
    lmeOut$coef$fixed[2],
    0.001)

omxCheckCloseEnough(mean(multilevelModel2Fit@output$estimate[(1:numSubjects)+(1*numSubjects)]), 
    lmeOut$coef$fixed[1],
    0.001)
