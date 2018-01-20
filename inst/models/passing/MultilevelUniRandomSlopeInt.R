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

require(OpenMx)
require(nlme)

# Multilevel Long Format Test
# Author: Steve Boker
# Date: Sun Nov 29 14:06:07 EST 2009


# This script is used to test the multilevel long format
# functionality using definition variables as indices.
totalOccasions <- 100
totalSubjects <- 10L
set.seed(42) # repeatibility
tID <- rep(1:totalSubjects, each=totalOccasions)
trueX <- rep(rnorm(totalOccasions, mean=0, sd=2), each=totalSubjects) +
    rnorm(totalOccasions*totalSubjects, mean=0, sd=.2)
trueB <- rep(rnorm(totalSubjects, mean=.8, sd=.3), each=totalOccasions)
tDataFrame <- data.frame(
    ID=tID, X=trueX, Y=trueB*trueX +
	rnorm(totalOccasions*totalSubjects,mean=0, sd=.1),trueB=trueB)
summary(tDataFrame)

manifestVars <- c("X", "Y")
numSubjects <- length(unique(tDataFrame$ID))

# Estimates the sum of the random and fixed effects
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
    mxFitFunctionML(),mxExpectationNormal(covariance="R", means="M"),
    mxData(tDataFrame, type="raw")
)

# ----------------------------------
# Fit the model and examine the summary results.

multilevelModel2Fit <- mxRun(multilevelModel2)

summary(multilevelModel2Fit)

lmeOut <- lme(Y~X, random= ~ X | ID, data=tDataFrame)

cbind(multilevelModel2Fit$output$estimate[1:numSubjects], 
      lmeOut$coef$random$ID[,2] + lmeOut$coef$fixed[2],
      trueB[seq(1,totalOccasions*(totalSubjects), by=totalOccasions)])

mean(multilevelModel2Fit$output$estimate[1:numSubjects])

est <- multilevelModel2Fit$output$estimate

omxCheckCloseEnough(mean(est[1:numSubjects]), 
    lmeOut$coef$fixed[2], 0.001)

omxCheckCloseEnough(mean(est[(1:numSubjects) + (1*numSubjects)]), 
    lmeOut$coef$fixed[1], 0.001)

# ----------------------------------
# An OpenMx equivalent to the mixed model

perID <- mxModel(
    "perID", type="RAM", latentVars=c('int', 'slope'),
    mxData(data.frame(ID=1L:totalSubjects), "raw", primaryKey="ID"),
    mxPath(c('int', 'slope'),c('int', 'slope'),'unique.pairs',
	   arrows=2,values=c(1,0,1)))

occa <- mxModel(
    "occa", type="RAM", perID, manifestVars="Y", latentVars="lX",
    mxData(tDataFrame, 'raw'),
    mxPath('Y', arrows=2, values=1),
    mxPath('one', 'Y'),
    mxPath('one', 'lX', labels='data.X', free=FALSE),
    mxPath('lX', 'Y'),
    mxPath('perID.int', 'Y', values=1, free=FALSE, joinKey='ID'),
    mxPath('perID.slope', 'Y', labels='data.X', free=FALSE, joinKey='ID'))

if (0) {
	require(lme4)
	lmer1 <- lmer(Y~X + (X | ID), data=tDataFrame, REML=FALSE)
	pt1 <- occa
	#pt1$perID$cholS$values[,] <- chol(VarCorr(lmer1)$ID)
	pt1$perID$S$values[,] <- VarCorr(lmer1)$ID
	pt1$A$values['Y', 'lX'] <- fixef(lmer1)['X']
	pt1$M$values[,'Y'] <- fixef(lmer1)['(Intercept)']
	pt1$S$values['Y', 'Y'] <- getME(lmer1, "sigma")^2

	pt1 <- mxRun(mxModel(pt1, mxComputeSequence(list(
	    mxComputeOnce('fitfunction', 'fit'),
	    mxComputeReportExpectation()))))

	omxCheckCloseEnough(logLik(pt1), logLik(lmer1), 1e-6)
}

occa <- mxRun(occa)
# a tad better than lme, same as lmer
omxCheckCloseEnough(occa$output$fit, -1725.954, 1e-2)
