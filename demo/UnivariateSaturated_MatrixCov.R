require(OpenMx)

#Simulate Data
set.seed(100)
x <- rnorm (1000, 0, 1)
testData <- as.matrix(x)
selVars <- c("X")
dimnames(testData) <- list(NULL, selVars)
summary(testData)
mean(testData)
var(testData)

#example 3: Saturated Model with Cov Matrices and Matrix-Style Input
univSatModel3 <- mxModel("univSat3",
    mxMatrix(
        type="Symm", 
        nrow=1, 
        ncol=1, 
        free=T, 
        values=1, 
        dimnames=list(selVars,selVars), 
        name="expCov"
    ),
    mxData(
        observed=var(testData), 
        type="cov", 
        numObs=1000 
    ),
    mxMLObjective(
        covariance="expCov"
    )
    )
univSatFit3 <- mxRun(univSatModel3)
EC3 <- mxEval(expCov, univSatFit3)
LL3 <- mxEval(objective, univSatFit3)
SL3 <- univSatFit3@output$SaturatedLikelihood
Chi3 <- LL3-SL3

#example 3m: Saturated Model with Cov Matrices & Means and Matrix-Style Input
univSatModel3m <- mxModel("univSat3m",
    mxMatrix(
        type="Symm", 
        nrow=1, 
        ncol=1, 
        free=T, 
        values=1, 
        dimnames=list(selVars,selVars), 
        name="expCov"
    ),
    mxMatrix(
        type="Full", 
        nrow=1, 
        ncol=1, 
        free=T, 
        values=0, 
        dimnames=list(selVars,NULL), 
        name="expMean"
    ),
    mxData(
        observed=var(testData), 
        type="cov", 
        numObs=1000,
        means=mean(testData) 
    ),
    mxMLObjective(
        covariance="expCov", 
        means="expMean"
    )
    )
univSatFit3m <- mxRun(univSatModel3m)
EM3m <- mxEval(expMean, univSatFit3m)
EC3m <- mxEval(expCov, univSatFit3m)
LL3m <- mxEval(objective, univSatFit3m);
SL3m <- univSatFit3m@output$SaturatedLikelihood
Chi3m <- LL3m-SL3m


#Mx answers hard-coded
#example Mx..1: Saturated Model with Cov Matrices
Mx.EC1 <-  1.062112
Mx.LL1 <- -1.474434e-17

#example Mx..1m: Saturated Model with Cov Matrices & Means
Mx.EM1m <- 0.01680509
Mx.EC1m <- 1.062112
Mx.LL1m <- -1.108815e-13


#Compare OpenMx results to Mx results (LL: likelihood; EC: expected covariance, EM: expected means)
##3:CovMat
omxCheckCloseEnough(Chi3,Mx.LL1,.001)
omxCheckCloseEnough(EC3,Mx.EC1,.001)
#3m:CovMPat 
omxCheckCloseEnough(Chi3m,Mx.LL1m,.001)
omxCheckCloseEnough(EC3m,Mx.EC1m,.001)
omxCheckCloseEnough(EM3m,Mx.EM1m,.001)
