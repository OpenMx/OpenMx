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

#example 3: Saturated Model with Cov Matrices and Matrix-Style Input
bivSatModel3 <- mxModel("bivSat3",
    mxMatrix(
        type="Symm", 
        nrow=2, 
        ncol=2, 
        free=T, 
        values=c(1,.5,1), 
        dimnames=list(selVars,selVars), 
        name="expCov"
    ),
    mxData(
        observed=cov(testData), 
        type="cov", 
        numObs=1000 
    ),
    mxMLObjective(
        covariance="expCov"
    )
    )
bivSatFit3 <- mxRun(bivSatModel3)
EC3 <- mxEval(expCov, bivSatFit3)
LL3 <- mxEval(objective,bivSatFit3)
SL3 <- summary(bivSatFit3)$SaturatedLikelihood
Chi3 <- LL3-SL3

#example 3m: Saturated Model with Cov Matrices & Means and Matrix-Style Input
bivSatModel3m <- mxModel("bivSat3m",
    mxMatrix(
        type="Symm", 
        nrow=2, 
        ncol=2, 
        free=T, 
        values=c(1,.5,1), 
        dimnames=list(selVars,selVars), 
        name="expCov"
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
        observed=cov(testData), 
        type="cov", 
        numObs=1000, 
        means=colMeans(testData) 
    ),
    mxMLObjective(
        covariance="expCov",
        means="expMean"
    )
    )
bivSatFit3m <- mxRun(bivSatModel3m)
EM3m <- mxEval(expMean, bivSatFit3m)
EC3m <- mxEval(expCov, bivSatFit3m)
LL3m <- mxEval(objective,bivSatFit3m)
SL3m <- summary(bivSatFit3m)$SaturatedLikelihood
Chi3m <- LL3m-SL3m


#Mx answers hard-coded
#example Mx..1: Saturated Model with Cov Matrices
Mx.EC1 <- matrix(c(1.0102951, 0.4818317, 0.4818317, 0.9945329),2,2)
Mx.LL1 <- -2.258885e-13

#example Mx..1m: Saturated Model with Cov Matrices & Means
Mx.EM1m <- matrix(c(0.03211648, -0.004883811),1,2)
Mx.EC1m <- matrix(c(1.0102951, 0.4818317, 0.4818317, 0.9945329),2,2)
Mx.LL1m <- -5.828112e-14


#Compare OpenMx results to Mx results (LL: likelihood; EC: expected covariance, EM: expected means)
#3:CovMat
omxCheckCloseEnough(Chi3,Mx.LL1,.001)
omxCheckCloseEnough(EC3,Mx.EC1,.001)
#3m:CovMPat 
omxCheckCloseEnough(Chi3m,Mx.LL1m,.001)
omxCheckCloseEnough(EC3m,Mx.EC1m,.001)
omxCheckCloseEnough(EM3m,Mx.EM1m,.001)
