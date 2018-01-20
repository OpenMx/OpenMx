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

#Simulate Data

set.seed(200)
rs=.5
xy <- mvtnorm::rmvnorm (1000, c(0,0), matrix(c(1,rs,rs,1),2,2))
testData <- xy
testData <- testData[, order(apply(testData, 2, var))[2:1]] #put the data columns in order from largest to smallest variance
# Note: Users do NOT have to re-order their data columns.  This is only to make data generation the same on different operating systems.
selVars <- c('X','Y')
dimnames(testData) <- list(NULL, selVars)
summary(testData)
cov(testData)

#Fit Saturated Model with RawData and Matrices Input using Cholesky Decomposition
bivCorModel <- mxModel("bivCor",
    mxMatrix(
        type="Full", 
        nrow=2, 
        ncol=2, 
        free=c(T,T,F,T), 
        values=c(1,.2,0,1), 
        name="Chol"
    ), 
    mxAlgebra(
        expression=Chol %*% t(Chol), 
        name="expCov"
    ), 
    mxMatrix(
        type="Full", 
        nrow=1, 
        ncol=2, 
        free=T, 
        values=c(0,0),
        name="expMean"
    ), 
    mxData(
        observed=testData, 
        type="raw"
    ), 
    mxFitFunctionML(),mxExpectationNormal(
        covariance="expCov", 
        means="expMean",
        dimnames=selVars
    )
    )

bivCorFit <- mxRun(bivCorModel)
EM <- mxEval(expMean, bivCorFit)
EC <- mxEval(expCov, bivCorFit)
LL <- mxEval(objective, bivCorFit);


#Test for Covariance=Zero
bivCorModelSub <-mxModel(bivCorModel,
    mxMatrix(
        type="Full", 
        nrow=2, 
        ncol=2, 
        free=c(T,F,F,T), 
        values=c(1,0,0,1), 
        name="Chol"
        )
    )
bivCorFitSub <- mxRun(bivCorModelSub)
EMs <- mxEval(expMean, bivCorFitSub)
ECs <- mxEval(expCov, bivCorFitSub)
LLs <- mxEval(objective, bivCorFitSub);
Chi= LLs-LL;
LRT= rbind(LL,LLs,Chi); LRT

#Compare OpenMx results to Mx results (LL: likelihood; EC: expected covariance, EM: expected means)
omxCheckCloseEnough(LL, 5407.036, .001)
omxCheckCloseEnough(c(EC), c(1.0656, 0.4752, 0.4752, 0.9292), .001)
omxCheckCloseEnough(c(EM), c(0.0582, 0.0063), .001)
