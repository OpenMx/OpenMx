# ===========
# = history =
# ===========
# 2017-04-14 04:35PM TBATES: update for mxFitFunctionMultigroup(c("MZ", "DZ"))
# CHECK stil FAILING at line 167
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

# j repeatedly gets set to a large value at subnp.cpp line 893.
#
# This is the same algorithmic problem with RSOLNP and CSOLNP we talked
# about some time ago. So, CSOLNP needs boundaries for some problems to
# get to the optimum. If you remove the boundary, it will get stuck in a
# loop.

if (mxOption(NULL,"Default optimizer") != 'CSOLNP') stop("SKIP")

# mxOption(NULL, "Default optimizer", "NPSOL")

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
		mxMatrix("Full", nrow=1, ncol=1, free=T, values=.6, lbound=.1, label="a", name="X"), 
		mxMatrix("Full", nrow=1, ncol=1, free=T, values=sqrt(.6), label="c", name="Y"),
		mxMatrix("Full", nrow=1, ncol=1, free=T, values=.6, lbound=.1, label="e", name="Z"),
		# Matrices A, C, and E compute variance components
		mxAlgebra(expression=X %*% t(X), name="A"), 
		mxAlgebra(expression=Y %*% t(Y), name="C"), 
		mxAlgebra(expression=Z %*% t(Z), name="E"),
		mxMatrix("Full", nrow=1, ncol=2, free=T, values= 20, lbound=0, ubound=30, label="mean", name="expMean"),
	    # Algebra for expected variance/covariance matrix in MZ
	    mxAlgebra(name="expCovMZ",
					rbind(cbind(A+C+E, A+C),
								cbind(A+C  , A+C+E))
		),
	    # Algebra for expected variance/covariance matrix in DZ
	    # note use of 0.5, converted to 1*1 matrix
	    mxAlgebra(name="expCovDZ", expression= 
				rbind(cbind(A+C+E     , 0.5%x%A+C),
							cbind(0.5%x%A+C , A+C+E))
			)
	),
	mxModel("MZ",
	    mxData(observed=mzfData, type="raw"), 
	    mxFitFunctionML(),
			mxExpectationNormal("common.expCovMZ", means="common.expMean", dimnames=selVars)
	),
	mxModel("DZ",
	    mxData(observed=dzfData, type="raw"), 
	    mxFitFunctionML(),
			mxExpectationNormal("common.expCovDZ", means="common.expMean", dimnames=selVars)
	),
	mxFitFunctionMultigroup(c("MZ", "DZ")),
  mxCI(c("common.A", "common.C[,]", "common.E[1,1]"))
)	
#Run ACE model
# -----------------------------------------------------------------------
twinACENoIntervals <- mxRun(twinACEModel, suppressWarnings = TRUE)
twinACENoIntervals <- mxRun(twinACEModel)
twinACEFit <- mxRun(twinACEModel, intervals=TRUE)

# summary(twinACEFit)
detail <- twinACEFit$compute$steps[['CI']]$output[['detail']]
print(detail)
omxCheckTrue(is.factor(detail[['side']]))
omxCheckEquals(levels(detail[['side']]), c('upper', 'lower'))

ci <- twinACEFit$output$confidenceIntervals
# cat(deparse(round(ci[,'ubound'],4)))
omxCheckCloseEnough(ci[-2,'lbound'], c(0.556, 0.1537), .01)
omxCheckCloseEnough(ci[,'ubound'], c(0.683, 0.052, 0.1956), .005)

iterateMxRun <- function(model, maxIterations) {
  model <- mxOption(model, "Optimality tolerance", 1e-6)
  return(iterateMxRunHelper(mxRun(mxModel(model, mxComputeGradientDescent(maxMajorIter=150L))),
				  maxIterations, 1))
}

iterateMxRunHelper <- function(model, maxIterations, iteration) {
	if (length(model$output) > 0 && model$output$status[[1]] == 0) {
		return(model)
	} else if (iteration < maxIterations) {
		return(iterateMxRunHelper(mxRun(mxModel(model, mxComputeGradientDescent(maxMajorIter=150L))),
					  maxIterations, iteration + 1))
	} else {
		return(model)
	}
}

twinACEIntervals <- twinACEFit
# twinACEIntervals$output <- list()

maxMisfit <- mxEval(objective, twinACEFit)[1,1] + 3.84

CIaupper <- mxModel(twinACEIntervals, name = 'A_CIupper',
    mxConstraint(MZ.objective + DZ.objective < maxMisfit),
		mxAlgebra(- common.A, name="upperCIa"), 
		mxFitFunctionAlgebra("upperCIa")
)

CIalower <- mxModel(twinACEIntervals, name = 'A_CIlower',
		mxConstraint(MZ.objective + DZ.objective < maxMisfit),
		mxAlgebra(common.A, name="lowerCIa"),
		mxFitFunctionAlgebra("lowerCIa")
)

runCIalower <- suppressWarnings(iterateMxRun(CIalower, 3))
runCIaupper <- suppressWarnings(iterateMxRun(CIaupper, 3))

CIcupper <- mxModel(twinACEIntervals, name = 'C_CIupper',
		mxConstraint(MZ.objective + DZ.objective < maxMisfit),
		mxAlgebra(- common.C,name="upperCIc"), 
		mxFitFunctionAlgebra("upperCIc")
)

runCIcupper <- suppressWarnings(iterateMxRun(CIcupper, 3))

CIeupper <- mxModel(twinACEIntervals, name = 'E_CIupper',
		mxConstraint(MZ.objective + DZ.objective < maxMisfit),
		mxAlgebra(- common.E,name="upperCIe"), 
		mxFitFunctionAlgebra("upperCIe")
)

CIelower <- mxModel(twinACEIntervals, name = 'E_CIlower',
		mxConstraint(MZ.objective + DZ.objective < maxMisfit),
		mxAlgebra(common.E,name="lowerCIe"),
		mxFitFunctionAlgebra("lowerCIe")
)

runCIelower <- suppressWarnings(iterateMxRun(CIelower, 3))
runCIeupper <- suppressWarnings(iterateMxRun(CIeupper, 3))

omxCheckCloseEnough(twinACEFit$output$confidenceIntervals[1, 'lbound'], mxEval(common.A, runCIalower), .01)
omxCheckCloseEnough(twinACEFit$output$confidenceIntervals[1, 'ubound'], mxEval(common.A, runCIaupper), .01)

# Can go either way
# omxCheckTrue(is.na(twinACEFit$output$confidenceIntervals[2, 'lbound']))

# Next check still failing under 2.7.9
omxCheckCloseEnough(twinACEFit$output$confidenceIntervals[2, 'ubound'], mxEval(common.C, runCIcupper), .001)

omxCheckCloseEnough(twinACEFit$output$confidenceIntervals[3, 'lbound'], mxEval(common.E, runCIelower), .005)
omxCheckCloseEnough(twinACEFit$output$confidenceIntervals[3, 'ubound'], mxEval(common.E, runCIeupper), .005)

twinACEParallel <- omxParallelCI(twinACENoIntervals)

if (0) {
  twinACEFit$output$confidenceIntervals - twinACEParallel$output$confidenceIntervals
}

mask <- !is.na(twinACEFit$output$confidenceIntervals)
omxCheckCloseEnough(twinACEFit$output$confidenceIntervals[mask],
	twinACEParallel$output$confidenceIntervals[mask], .001)
