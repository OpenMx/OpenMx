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


# -----------------------------------------------------------------------------
# Program: RowObjectiveFIMLBivariateSaturated.R
# Author: Hermine Maes 
# Date: 2009.08.01 
#
# Purpose: 
#      Demonstrate an mxFitFunctionRow function implementation of FIML
#
# ModelType: Saturated
# DataType: Simulated Continuous
# Field: None 
#
# RevisionHistory:
#	HermineMaes -- 2009.10.08 updated & reformatted
#	RossGore -- 2011.04.10 modified to implement FIML via mxFitFunctionRow
#	MikeHunter -- 2011.04.11 debugged Gore implementation above
#	MikeHunter -- 2011.05.03 Renamed from BivariateCorrelation.R to 
#                                        FIMLRowObjectiveBivariateCorrelation.R
#	MikeHunter -- 2011.05.05 modified to use omxSelect* functions
#	MikeHunter -- 2011.05.26 Adjusted spacing & comments for readability.
#   MikeHunter -- 2011.06.30 Adjusted model structure to match docs file in User's Guide and corrected formula error.
#      Hermine Maes -- 2014.11.04 piecewise specification
# -----------------------------------------------------------------------------

require(OpenMx)
require(MASS)
# Load Library
# -----------------------------------------------------------------------------

"%&%" <- OpenMx::"%&%"  # ensure we don't use the %&% from Matrix

set.seed(200)
rs <- .5
xy <- mvtnorm::rmvnorm (1000, c(0,0), matrix(c(1, rs, rs, 1), nrow=2, ncol=2))
testData <- as.data.frame(xy)
testData <- testData[, order(apply(testData, 2, var))[2:1]] #put the data columns in order from largest to smallest variance
# Note: Users do NOT have to re-order their data columns.  This is only to make data generation the same on different operating systems: to fix an inconsistency with the mvtnorm::rmvnorm function in the MASS package.
testVars <- c('X','Y')
names(testData) <- testVars
summary(testData)
cov(testData)
# Simulate Data: two standardized variables X & Y with correlation of .5


# -----------------------------------------------------------------------------

bivCorModelSpec <- mxModel(
    name="FIML Saturated Bivariate",
    mxData( observed=testData, type="raw" ),
    mxMatrix( type="Full", nrow=1, ncol=2, free=TRUE, values=c(0,0), name="expMean" ), 
    mxMatrix( type="Symm", nrow=2, ncol=2, name='expCov', values=c(.21, .2, .2, .21), free=TRUE )
)
diag(bivCorModelSpec$expCov$lbound) <- .01

bivCorFiltering <- mxModel(
    model=bivCorModelSpec,
    mxAlgebra( expression=omxSelectRowsAndCols(expCov, existenceVector), name="filteredExpCov" ),
    mxAlgebra( expression=omxSelectCols(expMean, existenceVector), name="filteredExpMean" ),
    mxAlgebra( expression=sum(existenceVector), name="numVar_i")
)

bivCorCalc <- mxModel(
    model=bivCorFiltering,
    mxAlgebra( expression=log(2*pi), name="log2pi" ),
    mxAlgebra( expression=log2pi %*% numVar_i + log(det(filteredExpCov)), name="firstHalfCalc" ),
    mxAlgebra( expression=(filteredDataRow - filteredExpMean) %&% solve(filteredExpCov), name="secondHalfCalc" )
)

bivCorRowObj <- mxModel(
    model=bivCorCalc,
    mxAlgebra( expression=(firstHalfCalc + secondHalfCalc), name="rowAlgebra" ),
    mxAlgebra( expression=sum(rowResults), name="reduceAlgebra" ),
    mxFitFunctionRow( rowAlgebra='rowAlgebra', reduceAlgebra='reduceAlgebra', dimnames=c('X','Y') )
)

bivCorTotal <- bivCorRowObj

# Fit Saturated Model with Raw Data and Matrix-style Input.
# Estimate with Full Information Maximum Likelihood (FIML).
# FIML is implemented here as an mxFitFunctionRow function for pedagogical reasons only.
# FIML for one row of data is
#   2*log(2*pi) + log(det(Cov)) + (Row - Mean) %*% solve(Cov) %*% t(Row - Mean)
#  where Cov is the filtered expected covariance matrix
#        Row is the filtered data row
#        Mean is the filtered expected means row vector
#        solve(*) is the inverse of *
#        t(*) is the transpose of *
#        det(*) is the determinant of *
# FIML for a whole data set is the sum of all the FIML rows.
# -----------------------------------------------------------------------------

bivCorFit <- mxRun(bivCorTotal)
EM <- mxEval(expMean, bivCorFit)
EC <- mxEval(expCov, bivCorFit)
LL <- mxEval(objective, bivCorFit)
# Run Model and Generate Output
# -----------------------------------------------------------------------------

omxCheckCloseEnough(LL, 5407.037, .001)
omxCheckCloseEnough(c(EC), c(1.066, 0.475, 0.475, 0.929), .001)
omxCheckCloseEnough(c(EM), c(0.058, 0.006), .001)
# Compare OpenMx Results to Mx Results 
# LL: likelihood; EC: expected covariance, EM: expected means
# -----------------------------------------------------------------------------
