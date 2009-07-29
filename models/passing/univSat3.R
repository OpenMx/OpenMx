require(OpenMx)

#Simulate Data
set.seed(100); x <- rnorm (1000, 0, 1)
testData <- as.matrix(x); selVars <- c("X"); dimnames(testData) <- list(NULL, selVars)
summary(testData); mean(testData); var(testData)


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


#Compare OpenMx results to Mx results (LL: likelihood; EC: expected covariance, EM: expected means)
#3:CovMat
omxCheckCloseEnough(Chi3,Mx.LL1,.001)
omxCheckCloseEnough(EC3,Mx.EC1,.001)
#3m:CovMPat 
omxCheckCloseEnough(Chi3m,Mx.LL1m,.001)
omxCheckCloseEnough(EC3m,Mx.EC1m,.001)
omxCheckCloseEnough(EM3m,Mx.EM1m,.001)
