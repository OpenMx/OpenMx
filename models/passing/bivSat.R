require(OpenMx)

#Simulate Data
require(MASS)
set.seed(200); rs=.5; xy <- mvrnorm (1000, c(0,0), matrix(c(1,rs,rs,1),2,2))
testData <- xy; selVars <- c('X','Y'); dimnames(testData) <- list(NULL, selVars)
summary(testData); cov(testData)

#example 1: Saturated Model with Cov Matrices and Paths input
bivSatModel1 <- mxModel("bivSat1",
	manifestVars= selVars,
	mxPath(from="X", to="Y", arrows=2, free=T, values=.2, lbound=.01, labels="covXY"),
	mxPath(from=c("X", "Y"), arrows=2, free=T, values=1, lbound=.01, labels=c("varX","varY")),
	mxData(cov(testData), type="cov", numObs=1000),
	type="RAM")
bivSatFit1 <- mxRun(bivSatModel1)
EC1 <- bivSatFit1[['S']]@values
LL1 <- mxEvaluate(objective,bivSatFit1)
SL1 <- bivSatFit1@output$other$Saturated
Chi1 <- LL1-SL1

#example 1m: Saturated Model with Cov Matrices & Means and Paths input
bivSatModel1m <- mxModel("bivSat1m",
	manifestVars= selVars,
	mxPath(from="X", to="Y", arrows=2, free=T, values=.2, lbound=.01, labels="covXY"),
	mxPath(from=c("X", "Y"), arrows=2, free=T, values=1, lbound=.01, labels=c("varX","varY")),
	mxPath(from="one", to=c("X", "Y"), arrows=1, free=T, values=.01, labels=c("meanX","meanY")),
	mxData(cov(testData), type="cov", numObs=1000, means=colMeans(testData)),
	type="RAM")
bivSatFit1m <- mxRun(bivSatModel1m)
EM1m <- bivSatFit1m[['M']]@values
EC1m <- bivSatFit1m[['S']]@values
LL1m <- mxEvaluate(objective,bivSatFit1m)
SL1m <- bivSatFit1m@output$other$Saturated
Chi1m <- LL1m-SL1m

#example 2: Saturated Model with Raw Data and Path input
bivSatModel2 <- mxModel("bivSat2",
	manifestVars= selVars,
	mxPath(from="X", to="Y", arrows=2, free=T, values=.2, lbound=.01, labels="covXY"),
	mxPath(from=c("X", "Y"), arrows=2, free=T, values=1, lbound=.01, labels=c("varX","varY")),
	mxData(testData, type="raw"),
	type="RAM")
bivSatFit2 <- mxRun(bivSatModel2)
EM2 <- bivSatFit2[['M']]@values
EC2 <- bivSatFit2[['S']]@values
LL2 <- mxEvaluate(objective,bivSatFit2)

#example 2s: Saturated Model with Raw Data and Path input built upon Cov/Means version
bivSatModel2s <- mxModel(bivSatModel1,
	mxData(testData, type="raw"),
	type="RAM")
bivSatFit2s <- mxRun(bivSatModel2s)
EM2s <- bivSatFit2s[['M']]@values
EC2s <- bivSatFit2s[['S']]@values
LL2s <- mxEvaluate(objective,bivSatFit2s)

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

#examples 4: Saturated Model with Raw Data and Matrices input
bivSatModel4 <- mxModel("bivSat4",
 	mxMatrix("Symm", nrow=2, ncol=2, free=T, values=c(1,.5,1), name="expCov",dimnames=list(selVars,selVars)),
 	mxMatrix("Full", nrow=1, ncol=2, free=T, values=c(0,0), name="expMean", dimnames=list(NULL, selVars)),
 	mxData(testData, type="raw"),
 	mxFIMLObjective("expCov", "expMean"))
bivSatFit4 <- mxRun(bivSatModel4)
EC4 <- bivSatFit4[['expCov']]@values
EM4 <- bivSatFit4[['expMean']]@values
LL4 <- mxEvaluate(objective,bivSatFit4)

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

#examples 6: Saturated Model with RawData and Matrices Input
bivSatModel6 <- mxModel("bivSat6",
 	mxMatrix("Full", nrow=2, ncol=2, free=c(T,T,F,T), values=c(1,.2,0,1), name="Chol", dimnames=list(selVars, selVars)),
	mxAlgebra(Chol %*% t(Chol), name="expCov", dimnames=list(selVars, selVars)),
 	mxMatrix("Full", nrow=1, ncol=2, free=T, values=c(0,0), name="expMean", dimnames=list(NULL, selVars)),
 	mxData(testData, type="raw"),
 	mxFIMLObjective("expCov", "expMean"))
bivSatFit6 <- mxRun(bivSatModel6)
EM6 <- bivSatFit6[['expMean']]@values
#EC6 <- bivSatFit6[['expCov']]@values
EC6 <- mxEvaluate(Chol %*% t(Chol),bivSatFit6)
LL6 <- mxEvaluate(objective,bivSatFit6)



original.directory <- getwd()
setwd('temp-files')

cat("*\n",file='star')
write.table(var(testData), file="cov", row.names=FALSE, quote=FALSE, col.names=FALSE)

if (.Platform$OS.type == "windows") {
	system("cmd /c cat star cov > testData.cov")
} else {
	system("cat star cov > testData.cov")
}

write.table(colMeans(testData), file="mea", row.names=FALSE, quote=FALSE, col.names=FALSE)

if (.Platform$OS.type == "windows") {
	system("cmd /c cat star mea > testData.mea")
} else {
	system("cat star mea > testData.mea")
}

testDataDF<-as.data.frame(testData)
write.table(testDataDF, file="testData.rec", row.names=FALSE, na=".", quote=FALSE, col.names=FALSE)

setwd(original.directory)

#example Mx..1: Saturated Model with Cov Matrices
mymatrices1 <- omxOriginalMx("mx-scripts/bivSatR1.mx", "temp-files")
#attach(mymatrices1) #matrixName groupNumber . jobNumber
Mx.EC1 <- mymatrices1$X3.1
Mx.LL1 <- mymatrices1$F3.1;

#example Mx..1m: Saturated Model with Cov Matrices & Means
mymatrices1m <- omxOriginalMx("mx-scripts/bivSatR1m.mx", "temp-files")
Mx.EM1m <- mymatrices1m$M3.1
Mx.EC1m <- mymatrices1m$X3.1
Mx.LL1m <- mymatrices1m$F3.1;

#example Mx..2: Saturated Model with Raw Data
mymatrices2 <- omxOriginalMx("mx-scripts/bivSatR2.mx", "temp-files")
Mx.EM2 <- mymatrices2$M3.1
Mx.EC2 <- mymatrices2$X3.1
Mx.LL2 <- mymatrices2$F3.1;

#OpenMx summary
cov <- rbind(cbind(EC1,EC1m,EC2),cbind(EC3,EC3m,EC4))
mean <- rbind(cbind(EM1m, EM2),cbind(EM3m,EM4))
like <- rbind(cbind(LL1,LL1m,LL2),cbind(LL3,LL3m,LL4))
cov; mean; like

#old Mx summary
Mx.cov <- cbind(Mx.EC1,Mx.EC1m,Mx.EC2)
Mx.mean <- cbind(Mx.EM1m,Mx.EM2)
Mx.like <- cbind(Mx.LL1,Mx.LL1m,Mx.LL2)
Mx.cov; Mx.mean; Mx.like


#Compare OpenMx results to Mx results (LL: likelihood; EC: expected covariance, EM: expected means)
#1:CovPat
omxCheckCloseEnough(Chi1,Mx.LL1,.001)
omxCheckCloseEnough(EC1,Mx.EC1,.001)
#1m:CovMPat 
omxCheckCloseEnough(Chi1m,Mx.LL1m,.001)
omxCheckCloseEnough(EC1m,Mx.EC1m,.001)
omxCheckCloseEnough(EM1m,Mx.EM1m,.001)
#2:RawPat 
omxCheckCloseEnough(LL2,Mx.LL2,.001)
omxCheckCloseEnough(EC2,Mx.EC2,.001)
omxCheckCloseEnough(EM2,Mx.EM2,.001)
#2:RawSPat
omxCheckCloseEnough(LL2s,Mx.LL2,.001)
omxCheckCloseEnough(EC2s,Mx.EC2,.001)
omxCheckCloseEnough(EM2s,Mx.EM2,.001)
#3:CovMat
omxCheckCloseEnough(Chi3,Mx.LL1,.001)
omxCheckCloseEnough(EC3,Mx.EC1,.001)
#3m:CovMPat 
omxCheckCloseEnough(Chi3m,Mx.LL1m,.001)
omxCheckCloseEnough(EC3m,Mx.EC1m,.001)
omxCheckCloseEnough(EM3m,Mx.EM1m,.001)
#4:RawMat
omxCheckCloseEnough(LL4,Mx.LL2,.001)
omxCheckCloseEnough(EC4,Mx.EC2,.001)
omxCheckCloseEnough(EM4,Mx.EM2,.001)

#5:CovMat Cholesky
omxCheckCloseEnough(Chi5,Mx.LL1,.001)
omxCheckCloseEnough(EC5,Mx.EC1,.001)
#5m:CovMPat Cholesky
omxCheckCloseEnough(Chi5m,Mx.LL1m,.001)
omxCheckCloseEnough(EC5m,Mx.EC1m,.001)
omxCheckCloseEnough(EM5m,Mx.EM1m,.001)
#6:RawMat Cholesky
omxCheckCloseEnough(LL6,Mx.LL2,.001)
omxCheckCloseEnough(EC6,Mx.EC2,.001)
omxCheckCloseEnough(EM6,Mx.EM2,.001)

