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

#example 2: Saturated Model with Raw Data and Path input
bivSatModel2 <- mxModel("bivSat2",
    manifestVars= selVars,
    mxPath(
        from=c("X", "Y"), 
        arrows=2, 
        free=T, 
        values=1, 
        lbound=.01, 
        labels=c("varX","varY")
    ),
    mxPath(
        from="X", 
        to="Y", 
        arrows=2, 
        free=T, 
        values=.2, 
        lbound=.01, 
        labels="covXY"
    ),
    mxData(
        observed=testData, 
        type="raw", 
    ),
    type="RAM"
    )
bivSatFit2 <- mxRun(bivSatModel2)
EM2 <- mxEvaluate(M, bivSatFit2)
EC2 <- mxEvaluate(S, bivSatFit2)
LL2 <- mxEvaluate(objective, bivSatFit2)
SL2 <- summary(bivSatFit2)$SaturatedLikelihood
Chi2 <- LL2-SL2


#Mx answers hard-coded
#example Mx..2: Saturated Model with Raw Data
Mx.EM2 <- matrix(c(0.03211188, -0.004889211),1,2)
Mx.EC2 <- matrix(c(1.0092891, 0.4813504, 0.4813504, 0.9935366),2,2)
Mx.LL2 <- 5415.772


#Compare OpenMx results to Mx results (LL: likelihood; EC: expected covariance, EM: expected means)
#2:RawPat 
omxCheckCloseEnough(LL2,Mx.LL2,.001)
omxCheckCloseEnough(EC2,Mx.EC2,.001)
omxCheckCloseEnough(EM2,Mx.EM2,.001)
