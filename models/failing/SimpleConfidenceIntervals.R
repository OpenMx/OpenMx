#
#   Copyright 2007-2010 The OpenMx Project
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
			values=.6,  
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
	    mxFIMLObjective(
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
	    mxFIMLObjective(
			covariance="common.expCovDZ", 
			means="common.expMean", 
			dimnames=selVars
		)
	),
	mxAlgebra(
		expression=MZ.objective + DZ.objective, 
		name="twin"
	), 
	mxAlgebraObjective("twin"),
    omxInterval(c("common.A[1,1]", "common.C[1,1]", "common.E[1,1]"), 3.84, 3.84)
)	
#Run ACE model
# -----------------------------------------------------------------------
twinACEFit <- mxRun(twinACEModel, intervals=TRUE)

CIaupper <- mxModel(twinACEFit, name = 'A_CIupper',
		mxMatrix("Full", values=mxEval(objective, twinACEFit), name="oldfit"), 
		mxAlgebra(((oldfit + 3.84) - (MZ.objective + DZ.objective))^2 - common.A,name="upperCIa"), 
		mxAlgebraObjective("upperCIa"))

CIalower <- mxModel(twinACEFit, name = 'A_CIlower',
		mxMatrix("Full", values=mxEval(objective, twinACEFit), name="oldfit"),
		mxAlgebra(((oldfit + 3.84) - (MZ.objective + DZ.objective))^2 + common.A,name="lowerCIa"),
		mxAlgebraObjective("lowerCIa"))

runCIalower <- mxRun(CIalower, intervals=FALSE, suppressWarnings = TRUE)
runCIaupper <- mxRun(CIaupper, intervals=FALSE, suppressWarnings = TRUE)

CIcupper <- mxModel(twinACEFit, name = 'C_CIupper',
		mxMatrix("Full", values=mxEval(objective, twinACEFit), name="oldfit"), 
		mxAlgebra(((oldfit + 3.84) - (MZ.objective + DZ.objective))^2 - common.C,name="upperCIc"), 
		mxAlgebraObjective("upperCIc"))

CIclower <- mxModel(twinACEFit, name = 'C_CIlower',
		mxMatrix("Full", values=mxEval(objective, twinACEFit), name="oldfit"),
		mxAlgebra(((oldfit + 3.84) - (MZ.objective + DZ.objective))^2 + common.C,name="lowerCIc"),
		mxAlgebraObjective("lowerCIc"))

runCIclower <- mxRun(CIclower, intervals=FALSE, suppressWarnings = TRUE)
runCIcupper <- mxRun(CIcupper, intervals=FALSE, suppressWarnings = TRUE)

CIeupper <- mxModel(twinACEFit, name = 'E_CIupper',
		mxMatrix("Full", values=mxEval(objective, twinACEFit), name="oldfit"), 
		mxAlgebra(((oldfit + 3.84) - (MZ.objective + DZ.objective))^2 - common.E,name="upperCIe"), 
		mxAlgebraObjective("upperCIe"))

CIelower <- mxModel(twinACEFit, name = 'E_CIlower',
		mxMatrix("Full", values=mxEval(objective, twinACEFit), name="oldfit"),
		mxAlgebra(((oldfit + 3.84) - (MZ.objective + DZ.objective))^2 + common.E,name="lowerCIe"),
		mxAlgebraObjective("lowerCIe"))

runCIelower <- mxRun(CIelower, intervals=FALSE, suppressWarnings = TRUE)
runCIeupper <- mxRun(CIeupper, intervals=FALSE, suppressWarnings = TRUE)


omxCheckCloseEnough(twinACEFit@output$confidenceIntervalCodes[1,1], runCIalower@output$status[1], .1)
omxCheckCloseEnough(twinACEFit@output$confidenceIntervalCodes[1,2], runCIaupper@output$status[1], .1)
omxCheckCloseEnough(twinACEFit@output$confidenceIntervals[1, 1], mxEval(common.A, runCIalower), .001)
omxCheckCloseEnough(twinACEFit@output$confidenceIntervals[1, 2], mxEval(common.A, runCIaupper), .001)

omxCheckCloseEnough(twinACEFit@output$confidenceIntervals[2, 1], mxEval(common.C, runCIclower), .001)
omxCheckCloseEnough(twinACEFit@output$confidenceIntervals[2, 2], mxEval(common.C, runCIcupper), .001)
omxCheckCloseEnough(twinACEFit@output$confidenceIntervalCodes[2,1], runCIclower@output$status[1], .1)
omxCheckCloseEnough(twinACEFit@output$confidenceIntervalCodes[2,2], runCIcupper@output$status[1], .1)


omxCheckCloseEnough(twinACEFit@output$confidenceIntervalCodes[3,1], runCIelower@output$status[1], .1)
omxCheckCloseEnough(twinACEFit@output$confidenceIntervalCodes[3,2], runCIeupper@output$status[1], .1)
omxCheckCloseEnough(twinACEFit@output$confidenceIntervals[3, 1], mxEval(common.E, runCIelower), .001)
omxCheckCloseEnough(twinACEFit@output$confidenceIntervals[3, 2], mxEval(common.E, runCIeupper), .001)

twinACEFit@output$confidenceIntervalCodes


