# -----------------------------------------------------------------------
# Program: UnivariateSaturated_PathRaw.R  
#  Author: Hermine Maes
#    Date: 08 01 2009 
#
# Univariate Saturated model to estimate means and variances
# Path style model input - Raw data input
#
# Revision History
#   Hermine Maes -- 10 08 2009 updated & reformatted
# -----------------------------------------------------------------------

require(OpenMx)

#Simulate Data
# -----------------------------------------------------------------------
set.seed(100)
x <- rnorm (1000, 0, 1)
testData <- as.matrix(x)
selVars <- c("X")
dimnames(testData) <- list(NULL, selVars)
summary(testData)
mean(testData)
var(testData)

#example 2: Saturated Model with Raw Data and Path-Style Input
# -----------------------------------------------------------------------
univSatModel2 <- mxModel("univSat2",
    manifestVars= selVars,
    mxPath(
        from=c("X"), 
        arrows=2, 
        free=T, 
        values=1, 
        lbound=.01, 
        labels="vX"
    ),
    mxData(
        observed=testData, 
        type="raw", 
    ),
    type="RAM"
)

univSatFit2 <- mxRun(univSatModel2)
EM2 <- mxEval(M, univSatFit2)
EC2 <- mxEval(S, univSatFit2)
LL2 <- mxEval(objective,univSatFit2);


#Mx answers hard-coded
# -----------------------------------------------------------------------
#example Mx..2: Saturated Model with Raw Data
Mx.EM2 <- 0.01680516
Mx.EC2 <- 1.061050
Mx.LL2 <- 2897.135


#Compare OpenMx results to Mx results
# -----------------------------------------------------------------------
# (LL: likelihood; EC: expected covariance, EM: expected means)
#2:RawPat 
omxCheckCloseEnough(LL2,Mx.LL2,.001)
omxCheckCloseEnough(EC2,Mx.EC2,.001)
omxCheckCloseEnough(EM2,Mx.EM2,.001)
