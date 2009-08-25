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

#examples 4: Saturated Model with Raw Data and Matrix-Style Input
univSatModel4 <- mxModel("univSat4",
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
        observed=testData, 
        type="raw", 
    ),
    mxFIMLObjective(
        covariance="expCov", 
        means="expMean"
    )
    )
univSatFit4 <- mxRun(univSatModel4)
EM4 <- mxEval(expMean, univSatFit4)
EC4 <- mxEval(expCov, univSatFit4)
LL4 <- mxEval(objective, univSatFit4);


#Mx answers hard-coded
#example Mx..1: Saturated Model with Raw Data
Mx.EM2 <- 0.01680516
Mx.EC2 <- 1.061050
Mx.LL2 <- 2897.135


#Compare OpenMx results to Mx results (LL: likelihood; EC: expected covariance, EM: expected means)
#4:RawMat
omxCheckCloseEnough(LL4,Mx.LL2,.001)
omxCheckCloseEnough(EC4,Mx.EC2,.001)
omxCheckCloseEnough(EM4,Mx.EM2,.001)

