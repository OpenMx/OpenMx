#
#   Copyright 2007-2015 The OpenMx Project
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

#Prepare Data
# -----------------------------------------------------------------------
data(twinData)
summary(twinData)
selVars <- c('bmi1','bmi2')
mzfData <- as.matrix(subset(twinData, zyg==1, c(bmi1,bmi2)))
dzfData <- as.matrix(subset(twinData, zyg==3, c(bmi1,bmi2)))
colMeans(mzfData,na.rm=TRUE)
colMeans(dzfData,na.rm=TRUE)
cov(mzfData,use="complete")
cov(dzfData,use="complete")

#Fit ACE Model with RawData and Matrices Input
# -----------------------------------------------------------------------
twinACEModel <- mxModel("twinACE",
	mxModel("common",
		# Matrices X, Y, and Z to store a, c, and e path coefficients
		mxMatrix(
			type="Full", 
			nrow=1, 
			ncol=1, 
			free=TRUE,  
			values=.6,  
			label="a", 
			name="X"
		), 
		mxMatrix(
			type="Full", 
			nrow=1, 
			ncol=1, 
			free=TRUE,  
			values=.6,  
			label="c", 
			name="Y"
		),
		mxMatrix(
			type="Full", 
			nrow=1, 
			ncol=1, 
			free=TRUE,  
			values=.6, lbound=1e-2,
			label="e", 
			name="Z"
		),
		# Matrices A, C, and E compute variance components
		mxAlgebra(
			expression=X %*% t(X), 
			name="A"
		), 
		mxAlgebra(
			expression=Y %*% t(Y), 
			name="C"
		), 
		mxAlgebra(
			expression=Z %*% t(Z), 
			name="E"
		),
		mxMatrix(
			type="Full", 
			nrow=1, 
			ncol=2, 
			free=TRUE, 
			values= 20,
			label="mean", 
			name="expMean"
		),
	    # Algebra for expected variance/covariance matrix in MZ
	    mxAlgebra(
			expression= rbind  (cbind(A+C+E , A+C),
								cbind(A+C   , A+C+E)), 
			name="expCovMZ"
		),
	    # Algebra for expected variance/covariance matrix in DZ
	    # note use of 0.5, converted to 1*1 matrix
	    mxAlgebra(
			expression= rbind  (cbind(A+C+E     , 0.5%x%A+C),
								cbind(0.5%x%A+C , A+C+E)), 
			name="expCovDZ"
		)
	),
	mxModel("MZ",
	    mxData(
			observed=mzfData, 
			type="raw"
		), 
	    mxFitFunctionML(),mxExpectationNormal(
			covariance="common.expCovMZ", 
			means="common.expMean", 
			dimnames=selVars
		)
	),
	mxModel("DZ",
	    mxData(
			observed=dzfData, 
			type="raw"
		), 
	    mxFitFunctionML(),mxExpectationNormal(
			covariance="common.expCovDZ", 
			means="common.expMean", 
			dimnames=selVars
		)
	),
	mxAlgebra(
		expression=MZ.objective + DZ.objective, 
		name="twin"
	), 
	mxFitFunctionAlgebra("twin"),
    mxCI(c("common.A", "common.C[,]", "common.E[1,1]"))
)	
#Run ACE model
# -----------------------------------------------------------------------
twinACENoIntervals <- mxRun(twinACEModel, suppressWarnings = TRUE)
twinACEFit <- mxRun(twinACEModel, intervals=TRUE, suppressWarnings = TRUE)

summary(twinACEFit)

iterateMxRun <- function(model, maxIterations) {
	return(iterateMxRunHelper(mxRun(model), maxIterations, 1))
}

iterateMxRunHelper <- function(model, maxIterations, iteration) {
	if (length(model$output) > 0 && model$output$status[[1]] == 0) {
		return(model)
	} else if (iteration < maxIterations) {
		return(iterateMxRunHelper(mxRun(model), maxIterations, iteration + 1))
	} else {
		return(model)
	}
}

twinACEIntervals <- twinACEFit
# twinACEIntervals$output <- list()

CIaupper <- mxModel(twinACEIntervals, name = 'A_CIupper',
		mxMatrix("Full", values=mxEval(objective, twinACEFit), name="oldfit"), 
		mxAlgebra(((oldfit + 3.84) - (MZ.objective + DZ.objective))^2 - common.A,name="upperCIa"), 
		mxFitFunctionAlgebra("upperCIa"))

CIalower <- mxModel(twinACEIntervals, name = 'A_CIlower',
		mxMatrix("Full", values=mxEval(objective, twinACEFit), name="oldfit"),
		mxAlgebra(((oldfit + 3.84) - (MZ.objective + DZ.objective))^2 + common.A,name="lowerCIa"),
		mxFitFunctionAlgebra("lowerCIa"))

runCIalower <- suppressWarnings(iterateMxRun(CIalower, 3))
runCIaupper <- suppressWarnings(iterateMxRun(CIaupper, 3))

CIcupper <- mxModel(twinACEIntervals, name = 'C_CIupper',
		mxMatrix("Full", values=mxEval(objective, twinACEFit), name="oldfit"), 
		mxAlgebra(((oldfit + 3.84) - (MZ.objective + DZ.objective))^2 - common.C,name="upperCIc"), 
		mxFitFunctionAlgebra("upperCIc"))

CIclower <- mxModel(twinACEIntervals, name = 'C_CIlower',
		mxMatrix("Full", values=mxEval(objective, twinACEFit), name="oldfit"),
		mxAlgebra(((oldfit + 3.84) - (MZ.objective + DZ.objective))^2 + common.C,name="lowerCIc"),
		mxFitFunctionAlgebra("lowerCIc"))

runCIclower <- suppressWarnings(iterateMxRun(CIclower, 3))
runCIcupper <- suppressWarnings(iterateMxRun(CIcupper, 3))

CIeupper <- mxModel(twinACEIntervals, name = 'E_CIupper',
		mxMatrix("Full", values=mxEval(objective, twinACEFit), name="oldfit"), 
		mxAlgebra(((oldfit + 3.84) - (MZ.objective + DZ.objective))^2 - common.E,name="upperCIe"), 
		mxFitFunctionAlgebra("upperCIe"))

CIelower <- mxModel(twinACEIntervals, name = 'E_CIlower',
		mxMatrix("Full", values=mxEval(objective, twinACEFit), name="oldfit"),
		mxAlgebra(((oldfit + 3.84) - (MZ.objective + DZ.objective))^2 + common.E,name="lowerCIe"),
		mxFitFunctionAlgebra("lowerCIe"))

runCIelower <- suppressWarnings(iterateMxRun(CIelower, 3))
runCIeupper <- suppressWarnings(iterateMxRun(CIeupper, 3))

omxCheckCloseEnough(twinACEFit$output$confidenceIntervals[1, 'lbound'], mxEval(common.A, runCIalower), .001)
omxCheckCloseEnough(twinACEFit$output$confidenceIntervals[1, 'ubound'], mxEval(common.A, runCIaupper), .001)

omxCheckCloseEnough(twinACEFit$output$confidenceIntervals[2, 'lbound'], mxEval(common.C, runCIclower), .001)
omxCheckCloseEnough(twinACEFit$output$confidenceIntervals[2, 'ubound'], mxEval(common.C, runCIcupper), .001)

omxCheckCloseEnough(twinACEFit$output$confidenceIntervals[3, 'lbound'], mxEval(common.E, runCIelower), .001)
omxCheckCloseEnough(twinACEFit$output$confidenceIntervals[3, 'ubound'], mxEval(common.E, runCIeupper), .001)

twinACEParallel <- omxParallelCI(twinACENoIntervals)

if (0) {
  twinACEFit$output$confidenceIntervals - twinACEParallel$output$confidenceIntervals
}

omxCheckCloseEnough(twinACEFit$output$confidenceIntervals, 
	twinACEParallel$output$confidenceIntervals, .001)
