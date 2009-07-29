require(OpenMx)

#Simulate Data
set.seed(100); x <- rnorm (1000, 0, 1)
testData <- as.matrix(x); selVars <- c("X"); dimnames(testData) <- list(NULL, selVars)
summary(testData); mean(testData); var(testData)

#example 1: Saturated Model with Cov Matrices and Paths input
univSatModel1 <- mxModel("univSat1",
	manifestVars= selVars,
	mxPath(from=c("X"), arrows=2, free=T, values=1, lbound=.01, labels="vX"),
	mxData(var(testData), type="cov", numObs=1000),
	type="RAM")
univSatFit1 <- mxRun(univSatModel1)
EC1 <- univSatFit1[['S']]@values
LL1 <- mxEvaluate(objective,univSatFit1)
SL1 <- univSatFit1@output$other$Saturated
Chi1 <- LL1-SL1

#example 1m: Saturated Model with Cov Matrices & Means and Paths input
univSatModel1m <- mxModel("univSat1m",
	manifestVars= selVars,
	mxPath(from=c("X"), arrows=2, free=T, values=.8, lbound=.01, labels="vX"),
	mxPath(from="one", to="X", arrows=1, free=T, values=0, labels="mX"),
	mxData(matrix(var(testData),1,1), type="cov", numObs=1000, means=mean(testData)),
	type="RAM")
univSatFit1m <- mxRun(univSatModel1m)
EM1m <- univSatFit1m[['M']]@values
EC1m <- univSatFit1m[['S']]@values
LL1m <- mxEvaluate(objective,univSatFit1m);
SL1m <- univSatFit1m@output$other$Saturated
Chi1m <- LL1m-SL1m

#example 2: Saturated Model with Raw Data and Path input
univSatModel2 <- mxModel("univSat2",
	manifestVars= selVars,
	mxPath(from=c("X"), arrows=2, free=T, values=.8, lbound=.01, labels="vX"),
	mxData(testData, type="raw"),
	type="RAM")
univSatFit2 <- mxRun(univSatModel2)
EM2 <- univSatFit2[['M']]@values
EC2 <- univSatFit2[['S']]@values
LL2 <- mxEvaluate(objective,univSatFit2);

#example 2s: Saturated Model with Raw Data and Path input built upon Cov/Means version
univSatModel2s <- mxModel(univSatModel1,
	mxData(testData, type="raw"),
	type="RAM")
univSatFit2s <- mxRun(univSatModel2s)
EM2s <- univSatFit2s[['M']]@values
EC2s <- univSatFit2s[['S']]@values
LL2s <- mxEvaluate(objective,univSatFit2s);

#example 3: Saturated Model with Cov Matrices and Matrices input
univSatModel3 <- mxModel("univSat3",
 	mxMatrix("Symm", nrow=1, ncol=1, free=T, values=1, name="expCov", dimnames=list(selVars,selVars)),
 	mxData(var(testData), type="cov", numObs=1000),
 	mxMLObjective("expCov"))
univSatFit3 <- mxRun(univSatModel3)
EC3 <- univSatFit3[['expCov']]@values
LL3 <- mxEvaluate(objective,univSatFit3);
SL3 <- univSatFit3@output$other$Saturated
Chi3 <- LL3-SL3

#example 3m: Saturated Model with Cov Matrices & Means and Matrices input
univSatModel3m <- mxModel("univSat3m",
 	mxMatrix("Symm", nrow=1, ncol=1, free=T, values=1, name="expCov", dimnames=list(selVars,selVars)),
 	mxMatrix("Full", nrow=1, ncol=1, free=T, values=0, name="expMean", dimnames=list(NULL, selVars)),
 	mxData(var(testData), type="cov", numObs=1000, means=mean(testData)),
 	mxMLObjective("expCov","expMean"))
univSatFit3m <- mxRun(univSatModel3m)
EM3m <- univSatFit3m[['expMean']]@values
EC3m <- univSatFit3m[['expCov']]@values
LL3m <- mxEvaluate(objective,univSatFit3m);
SL3m <- univSatFit3m@output$other$Saturated
Chi3m <- LL3m-SL3m

#examples 4: Saturated Model with Raw Data and Matrices input
univSatModel4 <- mxModel("univSat4",
 	mxMatrix("Symm", nrow=1, ncol=1, free=T, values=1, name="expCov",dimnames=list(selVars,selVars)),
 	mxMatrix("Full", nrow=1, ncol=1, free=T, values=0, name="expMean", dimnames=list(NULL, selVars)),
 	mxData(testData, type="raw"),
 	mxFIMLObjective("expCov", "expMean"))
univSatFit4 <- mxRun(univSatModel4)
EC4 <- univSatFit4[['expCov']]@values
EM4 <- univSatFit4[['expMean']]@values
LL4 <- mxEvaluate(objective,univSatFit4);

original.directory <- getwd()
setwd('temp-files')

cat("*\n",file='star')
write.table(var(testData), file="cov", row.names=FALSE, quote=FALSE, col.names=FALSE)

if (.Platform$OS.type == "windows") {
	system("cmd /c cat star cov > testData.cov")
} else {
	system("cat star cov > testData.cov")
}

write.table(mean(testData), file="mea", row.names=FALSE, quote=FALSE, col.names=FALSE)

if (.Platform$OS.type == "windows") {
	system("cmd /c cat star mea > testData.mea")
} else {
	system("cat star mea > testData.mea")
}

testDataDF<-as.data.frame(testData)
write.table(testDataDF, file="testData.rec", row.names=FALSE, na=".", quote=FALSE, col.names=FALSE)

setwd(original.directory)

#example Mx..1: Saturated Model with Cov Matrices
mymatrices1 <- omxOriginalMx("mx-scripts/univSatR1.mx", "temp-files")
#attach(mymatrices1) #matrixName groupNumber . jobNumber
Mx.EC1 <- mymatrices1$X3.1
Mx.LL1 <- mymatrices1$F3.1;

#example Mx..1m: Saturated Model with Cov Matrices & Means
mymatrices1m <- omxOriginalMx("mx-scripts/univSatR1m.mx", "temp-files")
Mx.EM1m <- mymatrices1m$M3.1
Mx.EC1m <- mymatrices1m$X3.1
Mx.LL1m <- mymatrices1m$F3.1;

#example Mx..2: Saturated Model with Raw Data
mymatrices2 <- omxOriginalMx("mx-scripts/univSatR2.mx", "temp-files")
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

