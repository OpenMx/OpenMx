# Saturated Model with RawData and Matrix-Style Input
# Based on Example 6 ...?
# Can we embed a page reference or clickable file:// or http:// link here?


require(OpenMx)
#Simulate Data
require(MASS)
set.seed(200)
rs=.5 # set the correlation of the simulated data
testData <- mvrnorm (1000, c(0,0), matrix(c(1,rs,rs,1),2,2)) # make the simulated data
selVars <- c('X','Y')
dimnames(testData) <- list(NULL, selVars)
summary(testData)
R.cov = cov(testData)

bivSatModel6 <- mxModel("bivSat6",
    mxMatrix(type="Full", nrow=2, ncol=2, free=c(TRUE,TRUE,FALSE,TRUE), values=c(1,.2,0,1), name="Chol"),
		mxMatrix(type="Full", nrow=1, ncol=2, free=TRUE, values=c(0,0), dimnames=list(NULL,selVars), name="expMean"),
    mxAlgebra(Chol %*% t(Chol), dimnames=list(selVars,selVars), name="expCov"),
		mxData(observed=testData, type="raw"),
		mxFIMLObjective(covariance="expCov",means="expMean")
		)
bivSatFit6 <- mxRun(bivSatModel6)
EM <- mxEval(expMean, bivSatFit6)
EC <- mxEval(expCov, bivSatFit6)
LL <- mxEval(objective,bivSatFit6)

#Mx answers hard-coded
#example Mx..2: Saturated Model with Raw Data
Mx.EM <- matrix(c(0.03211188, -0.004889211),1,2)
Mx.EC <- matrix(c(1.0092891, 0.4813504, 0.4813504, 0.9935366),2,2)
Mx.LL <- 5415.772

# Compare OpenMx results to Mx results (LL: likelihood; EC: expected covariance, EM: expected means)
#6: RawMat Cholesky
omxCheckCloseEnough(LL,Mx.LL,.001)
omxCheckCloseEnough(EC,Mx.EC,.001)
omxCheckCloseEnough(EM,Mx.EM,.001)

# Print out the covariance we generated in R, and that from this model in OpenMx
R.cov
EC