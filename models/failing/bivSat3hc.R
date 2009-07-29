require(OpenMx)

#Simulate Data
require(MASS)
set.seed(200); rs=.5; xy <- mvrnorm (1000, c(0,0), matrix(c(1,rs,rs,1),2,2))
testData <- xy; selVars <- c('X','Y'); dimnames(testData) <- list(NULL, selVars)
summary(testData); cov(testData)


#example 3: Saturated Model with Cov Matrices and Matrices input
bivSatModel3 <- mxModel("bivSat3",
 	mxMatrix("Symm", nrow=2, ncol=2, free=T, values=c(1,.5,1), name="expCov", dimnames=list(selVars,selVars)),
 	mxData(cov(testData), type="cov", numObs=1000),
 	mxMLObjective("expCov"))
bivSatFit3 <- mxRun(bivSatModel3)
EC3 <- bivSatFit3[['expCov']]@values
LL3 <- mxEvaluate(objective,bivSatFit3)
SL3 <- bivSatFit3@output$other$Saturated
Chi3 <- LL3-SL3

#example 3m: Saturated Model with Cov Matrices & Means and Matrices input
bivSatModel3m <- mxModel("bivSat3m",
 	mxMatrix("Symm", nrow=2, ncol=2, free=T, values=c(1,.5,1), name="expCov", dimnames=list(selVars,selVars)),
 	mxMatrix("Full", nrow=1, ncol=2, free=T, values=c(0,0), name="expMean", dimnames=list(NULL, selVars)),
 	mxData(cov(testData), type="cov", numObs=1000, means=colMeans(testData)),
 	mxMLObjective("expCov","expMean"))
bivSatFit3m <- mxRun(bivSatModel3m)
EM3m <- bivSatFit3m[['expMean']]@values
EC3m <- bivSatFit3m[['expCov']]@values
LL3m <- mxEvaluate(objective,bivSatFit3m)
SL3m <- bivSatFit3m@output$other$Saturated
Chi3m <- LL3m-SL3m


#example 5: Saturated Model with Cov Matrices and Matrices Input
bivSatModel5 <- mxModel("bivSat5",
 	mxMatrix("Full", nrow=2, ncol=2, free=c(T,T,F,T), values=c(1,.2,0,1), name="Chol"),
	mxAlgebra(Chol %*% t(Chol), name="expCov", dimnames=list(selVars,selVars)),
 	mxData(cov(testData), type="cov", numObs=1000),
 	mxMLObjective("expCov"))
bivSatFit5 <- mxRun(bivSatModel5)
#not working: EC5 <- bivSatFit5[['expCov']]@values  ?? $algebras$bivSat5.expCov
EC5 <- mxEvaluate(Chol %*% t(Chol),bivSatFit5)
LL5 <- mxEvaluate(objective,bivSatFit5)
SL5 <- bivSatFit5@output$other$Saturated
Chi5 <- LL5-SL5

#example 5m: Saturated Model with Cov Matrices & Means and Matrices Input
bivSatModel5m <- mxModel("bivSat5m",
 	mxMatrix("Full", nrow=2, ncol=2, free=c(T,T,F,T), values=c(1,.2,0,1), name="Chol"),
	mxAlgebra(Chol %*% t(Chol), name="expCov", dimnames=list(selVars,selVars)),
 	mxMatrix("Full", nrow=1, ncol=2, free=T, values=c(0,0), name="expMean", dimnames=list(NULL, selVars)),
 	mxData(cov(testData), type="cov", numObs=1000, means=colMeans(testData)),
 	mxMLObjective("expCov","expMean"))
bivSatFit5m <- mxRun(bivSatModel5m)
EM5m <- bivSatFit5m[['expMean']]@values
#EC5m <- bivSatFit5m[['expCov']]@values
EC5m <- mxEvaluate(Chol %*% t(Chol),bivSatFit5m)
LL5m <- mxEvaluate(objective,bivSatFit5m);
SL5m <- bivSatFit5m@output$other$Saturated
Chi5m <- LL5m-SL5m

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

#5:CovMat Cholesky
omxCheckCloseEnough(Chi5,Mx.LL1,.001)
omxCheckCloseEnough(EC5,Mx.EC1,.001)
#5m:CovMPat Cholesky
omxCheckCloseEnough(Chi5m,Mx.LL1m,.001)
omxCheckCloseEnough(EC5m,Mx.EC1m,.001)
omxCheckCloseEnough(EM5m,Mx.EM1m,.001)
