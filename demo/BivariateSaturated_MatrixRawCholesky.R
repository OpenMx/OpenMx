require(OpenMx)

#Simulate Data
require(MASS)
set.seed(200)
rs=.5
xy <- mvrnorm (1000, c(0,0), matrix(c(1,rs,rs,1),2,2))
testData <- xy
selVars <- c('X','Y')
dimnames(testData) <- list(NULL, selVars)
summary(testData)
cov(testData)

#examples 6: Saturated Model with RawData and Matrix-Style Input
bivSatModel6 <- mxModel("bivSat6",
    mxMatrix(
        type="Full", 
        nrow=2, 
        ncol=2, 
        free=c(T,T,F,T), 
        values=c(1,.2,0,1), 
        name="Chol"
    ),
    mxAlgebra(
        Chol %*% t(Chol), 
        name="expCov", 
        dimnames=list(selVars,selVars)
    ),
    mxMatrix(
        type="Full", 
        nrow=1, 
        ncol=2, 
        free=T, 
        values=c(0,0), 
        dimnames=list(NULL,selVars), 
        name="expMean"
    ),
    mxData(
        observed=testData, 
        type="raw", 
    ),
    mxFIMLObjective(
        covariance="expCov",
        means="expMean"
    )
    )
bivSatFit6 <- mxRun(bivSatModel6)
EM6 <- mxEval(expMean, bivSatFit6)
EC6 <- mxEval(expCov, bivSatFit6)
LL6 <- mxEval(objective,bivSatFit6)


#Mx answers hard-coded
#example Mx..2: Saturated Model with Raw Data
Mx.EM2 <- matrix(c(0.03211188, -0.004889211),1,2)
Mx.EC2 <- matrix(c(1.0092891, 0.4813504, 0.4813504, 0.9935366),2,2)
Mx.LL2 <- 5415.772


#Compare OpenMx results to Mx results (LL: likelihood; EC: expected covariance, EM: expected means)
#6:RawMat Cholesky
omxCheckCloseEnough(LL6,Mx.LL2,.001)
omxCheckCloseEnough(EC6,Mx.EC2,.001)
omxCheckCloseEnough(EM6,Mx.EM2,.001)

