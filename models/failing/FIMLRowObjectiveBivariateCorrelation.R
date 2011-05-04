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


# -----------------------------------------------------------------------
# Program: FIMLRowObjectiveBivariateCorrelation.R
#  Author: Hermine Maes -- See Revision History
#    Date: 08 01 2009 
#
# Optimization Example in OpenMx: Testing significance of correlation
#
# Revision History
#   Hermine Maes -- 10 08 2009 updated & reformatted
#   Ross Gore -- 2011.04.10 modified to implement FIML via mxRowObjective
#   Mike Hunter -- 2011.04.11 debugged Gore implementation above
#   Mike Hunter -- 2011.05.03 Renamed from BivariateCorrelation.R to FIMLRowObjectiveBivariateCorrelation.R
# -----------------------------------------------------------------------

require(OpenMx)

# Simulate Data: two standardized variables X & Y with correlation of .5
# -----------------------------------------------------------------------
require(MASS)
set.seed(200)
rs=.5
xy <- mvrnorm (1000, c(0,0), matrix(c(1,rs,rs,1),2,2))
testData <- xy
selVars <- c('X','Y')
dimnames(testData) <- list(NULL, selVars)
summary(testData)
cov(testData)

# Make some mxAlgebras to test mxRowObjective
# -----------------------------------------------------------------------


# Fit Saturated Model with Raw Data and Matrix-style Input
# -----------------------------------------------------------------------
# I (Michael Hunter) have edited this model to not use
# 	omxSelect* functions.  In this condition it works wonderfully,
# 	but cannot handle missing data.
# As of Mon Apr 11 19:52:21 EDT 2011 I believe the problem is with
# 	the mxRowObjective populating the existenceVector.
bivCorModel <- mxModel("bivCor",
    mxMatrix(
        type="Full", 
        nrow=1, 
        ncol=2, 
        free=TRUE, 
        values=c(0,0), 
        name="expMean"
    ), 
    mxMatrix(
        type="Lower", 
        nrow=2, 
        ncol=2, 
        free=TRUE,
        values=.2, 
        name="Chol"
    ), 
    mxAlgebra(
        expression=Chol %*% t(Chol), 
        name="expCov", 
    ),
	mxMatrix("Full", 1, 1, values = log(2*pi), name = "log2pi"),
 	#mxAlgebra(
	#	expression=omxSelectRowsAndCols(expCov, existenceVector),
	#	name="filteredExpCov",
	#),
	#mxAlgebra(
	#	expression=omxSelectCols(expMean, existenceVector), #existenceVector
	#	name="filteredExpMean",
	#),
	mxAlgebra(
		expression=log2pi %*% 2 + log(det(expCov)), #filteredExpCov
		name ="firstHalfCalc",
	),
	mxAlgebra(
		expression=(filteredDataRow - expMean) %&% solve(expCov), #filteredExpCov
		name = "secondHalfCalc",
	),
	mxAlgebra(
		expression=(firstHalfCalc+secondHalfCalc),
		name="rowAlgebra",
	),
	mxAlgebra(
		expression=sum(rowResults),
		name = "reduceAlgebra",
	),
	mxRowObjective(
		rowAlgebra='rowAlgebra',
		reduceAlgebra='reduceAlgebra',
		dimnames=c('X','Y'),
	),
    mxData(
        observed=testData, 
        type="raw",
    )
)	


# Run Model and Generate Output
# -----------------------------------------------------------------------
bivCorFit <- mxRun(bivCorModel)
EM <- mxEval(expMean, bivCorFit)
EC <- mxEval(expCov, bivCorFit)
LL <- mxEval(objective, bivCorFit)


# Specify SubModel testing Covariance=Zero
# -----------------------------------------------------------------------
bivCorModelSub <-mxModel(bivCorModel, name = "bivCorSub",
    mxMatrix(
        type="Diag", 
        nrow=2, 
        ncol=2,
		values=c(1,1), 
        free=TRUE, 
        name="Chol"
    )
)

# Run Model and Generate Output
# -----------------------------------------------------------------------
bivCorFitSub <- mxRun(bivCorModelSub)
EMs <- mxEval(expMean, bivCorFitSub)
ECs <- mxEval(expCov, bivCorFitSub)
LLs <- mxEval(objective, bivCorFitSub)
Chi= LLs-LL;
LRT= rbind(LL,LLs,Chi); LRT


# Mx Answers of Saturated Model Hard-coded
# -----------------------------------------------------------------------
Mx.EM <- matrix(c(0.03211656, -0.004883885), 1, 2)
Mx.EC <- matrix(c(1.0092853, 0.4813504, 0.4813504, 0.9935390), 2, 2)
Mx.LL <- 5415.772

# Compare OpenMx Results to Mx Results 
# -----------------------------------------------------------------------
#LL: likelihood; EC: expected covariance, EM: expected means
omxCheckCloseEnough(LL,Mx.LL,.001)
omxCheckCloseEnough(EC,Mx.EC,.001)
omxCheckCloseEnough(EM,Mx.EM,.001)
