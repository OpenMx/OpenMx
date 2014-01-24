#
#   Copyright 2007-2009 The OpenMx Project
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

# Multilevel Long Format Test
# Author: Steve Boker
# Date: Sun Nov 29 14:06:07 EST 2009


# This script is used to test the multilevel long format functionality using definition variables
#    as indices.

set.seed(42) # repeatibility
tID <- rep(1:10, each=10)
trueX <- rep(rnorm(10, mean=0, sd=2), each=10) + rnorm(100, mean=0, sd=.2)
trueB <- rep(rnorm(10, mean=0, sd=.5), each=10)
tDataFrame <- data.frame(ID=tID, X=trueX, Y=trueB*trueX + rnorm(100, mean=0, sd=.2))
summary(tDataFrame)

manifestVars <- c("X", "Y")
numSubjects <- length(unique(tDataFrame$ID))

RandVals <- matrix(c(.2, 0, 0, .8, .8), numSubjects, 5, byrow=TRUE)

multilevelModel2 <- mxModel("Multilevel_2",
    mxMatrix("Full", 
        values=RandVals,
        free=TRUE, 
        name="Rand"
    ),
    mxMatrix("Full", 2, 2, 
        labels=c(NA,  NA,
                "Rand[data.ID,1]", NA), 
        free=FALSE, 
        name="A", 
        byrow=TRUE
    ),
    mxMatrix("Symm", 2, 2,
        free=FALSE, 
        labels=c("Rand[data.ID,4]",
                 NA, "Rand[data.ID,5]"), 
        name="S", 
        byrow=TRUE
    ),
	mxMatrix("Full", 2, 2, 
	         values=c(1,0,
		              0,1),
		     free=F, byrow=TRUE, name="F"),
    mxMatrix("Iden", 2, name="I"),
    mxAlgebra(F %*% solve(I-A) %*% S %*% t(solve(I-A)) %*% t(F), 
        name="R", 
        dimnames = list(manifestVars, manifestVars)
    ),
    mxMatrix("Full", nrow=1, ncol=length(manifestVars),
        values=0,
        free=FALSE,
        labels=c("Rand[data.ID,2]", "Rand[data.ID,3]"),
        dimnames=list(NULL, manifestVars),
        name="M"
    ),
    mxFIMLObjective("R", "M"),
    mxData(tDataFrame, 
        type="raw"
    )
)

# ----------------------------------
# Fit the  model and examine the summary results.

multilevelModel2Fit <- mxRun(multilevelModel2)

summary(multilevelModel2Fit)

