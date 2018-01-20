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


#  OpenMx Out of order thresholds test
require(OpenMx)

set.seed(1234)
# Set up simulation parameters 
nVariables  <- 3
varNames    <- paste0("var",1:nVariables)
nFactors    <- 1
nThresholds <- 3
nSubjects   <- 500

# Simulate multivariate normal data and chop into ordinal
loadings  <- matrix(.7,nrow=nVariables,ncol=nFactors)
residuals <- 1 - (loadings * loadings)
sigma     <- loadings %*% t(loadings) + vec2diag(residuals)
mu        <- matrix(0,nrow=nVariables,ncol=1)

continuousData <- mvtnorm::rmvnorm(n=nSubjects, mu, sigma)
quants <-quantile(continuousData[,1],  probs = c((1:nThresholds)/(nThresholds+1)))
quants[2] <- .65
ordinalData <- matrix(0, nrow = nSubjects, ncol = nVariables)
for(i in 1:nVariables) {
	ordinalData[,i] <- cut(as.vector(continuousData[,i]), c(-Inf, quants, Inf))
}
ordinalData <- mxFactor(as.data.frame(ordinalData),levels=c(1:(nThresholds+1)))
names(ordinalData) <- varNames
# table(list(ordinalData[,1],ordinalData[,2]))

m1 <- mxModel("m1",
	mxMatrix(name = "vectorofOnes", "Unit", nVariables, 1),
	mxMatrix(name = "L"           , "Full", nVariables, nFactors, values=0.2, free=TRUE, lbound=-.99, ubound=.99),
	mxMatrix(name = "M"           , "Zero", 1, nVariables),
	mxAlgebra(vectorofOnes - (diag2vec(L %*% t(L))) , name = "E"),
	mxAlgebra(L %*% t(L) + vec2diag(E), name="impliedCovs"),
	mxMatrix(name="thresholdDeviations", "Full", nrow=nThresholds, ncol=nVariables,
            values = c(.2,.201,.3),
            free = TRUE,
            lbound = rep(c(-Inf,rep(.01,(nThresholds-1))) , nVariables),
            dimnames = list(c(), varNames)
	),
    # mxMatrix("Lower",nThresholds,nThresholds,values=1,free=F,name="unitLower"),
    mxAlgebra(thresholdDeviations, name = "thresholdMatrix"),
	mxExpectationNormal("impliedCovs", means = "M", dimnames = varNames, thresholds = "thresholdMatrix"),
	mxFitFunctionML(),
	mxData(ordinalData, type = 'raw')
)
# set up checkpointing to observe the threshold locations
m1 <- mxOption(m1,'Checkpoint Units', 'evaluations')
m1 <- mxOption(m1,'Checkpoint Count', 1)
mxRun(m1, checkpoint=TRUE, unsafe= TRUE)
summary(m1)
