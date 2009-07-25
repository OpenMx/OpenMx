require(OpenMx)
setwd('/Users/hermine/Applications/bin/OpenMx/trunk/demo')
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
EC1 <- univSatFit1[['S']]@values; LL1 <- mxEvaluate(objective,univSatFit1);

#example 1m: Saturated Model with Cov Matrices, Means and Paths input
#univSatModel1m <- mxModel("univSat1m",
#	manifestVars= selVars,
#	mxPath(from=c("X"), arrows=2, free=T, values=.8, lbound=.01, labels="vX"),
#	mxPath(from="one" to "X", arrows=1, free=T, values=0, labels="mX"),
#	mxData(matrix(var(testData),1,1), type="cov", means=mean(testData), numObs=1000),
#	type="RAM")
#univSatFit1m <- mxRun(univSatModel1m)
#EC1m <- univSatFit1m[['S']]@values; LL1m <- mxEvaluate(objective,univSatFit1m);

#example 2: Saturated Model with Raw Data and Path input
univSatModel2 <- mxModel("univSat2",
	manifestVars= selVars,
	mxPath(from=c("X"), arrows=2, free=T, values=.8, lbound=.01, labels="vX"),
	mxData(testData, type="raw"),
	type="RAM")
univSatFit2 <- mxRun(univSatModel2)
EM2 <- univSatFit2[['M']]@values; EC2 <- univSatFit2[['S']]@values; LL2 <- mxEvaluate(objective,univSatFit2);

#univSatModel2s <- mxModel(univSatModel1,
#	mxData(as.matrix(testData), type="raw"),
#	type="RAM")
#univSatFit2s <- mxRun(univSatModel2s)
#EC2s <- univSatFit2s[['S']]@values; LL2s <- mxEvaluate(objective,univSatFit2s);

#example 3: Saturated Model with Cov Matrices and Matrices input
univSatModel3 <- mxModel("univSat3",
 	mxMatrix("Symm", nrow=1, ncol=1, free=T, values=1, name="expCov", dimnames=list(selVars,selVars)),
 	mxData(var(testData), type="cov", numObs=1000),
 	mxMLObjective("expCov"))
univSatFit3 <- mxRun(univSatModel3)
EC3 <- univSatFit3[['expCov']]@values; LL3 <- mxEvaluate(objective,univSatFit3);

#examples 4: Saturated Model with Raw Data and Matrices input
univSatModel4 <- mxModel("univSat4",
 	mxMatrix("Symm", nrow=1, ncol=1, free=T, values=1, name="expCov",dimnames=list(selVars,selVars)),
 	mxMatrix("Full", nrow=1, ncol=1, free=T, values=0, name="expMean", dimnames=list(NULL, selVars)),
 	mxData(testData, type="raw"),
 	mxFIMLObjective("expCov", "expMean"))
univSatFit4 <- mxRun(univSatModel4)
EC4 <- univSatFit4[['expCov']]@values; EM4 <- univSatFit4[['expMean']]@values; LL4 <- mxEvaluate(objective,univSatFit4);

setwd("/Users/hermine/Applications/bin/OpenMx/trunk/demo/MxR")
cat("*\n",file='star')
write.table(var(testData),file="cov",row.names=F,quote=F,col.names=F); system("cat star cov > testData.cov")
write.table(mean(testData),file="mea",row.names=F,quote=F,col.names=F); system("cat star mea > testData.mea")
testDataDF<-as.data.frame(testData)
write.table(testDataDF,file="testData.rec",row.names=F,na=".",quote=F,col.names=F)

source("runmx.R")
#Mx file includes 3 Saturated Model jobs: 1: Cov Matrices, 2: Cov Matrices + Means, 3: Raw Data 
mymatrices1 <- runmx(mx.script="univSatR1.mx",mx.location="/usr/local/bin/mxt.169")
attach(mymatrices1) #matrixName groupNumber . jobNumber
#example Mx..1: Saturated Model with Cov Matrices
Mx.EC1 <-X3.1; Mx.LL1 <- F3.1;

mymatrices2 <- runmx(mx.script="univSatR1m.mx",mx.location="/usr/local/bin/mxt.169")
attach(mymatrices2) #matrixName groupNumber . jobNumber
#example Mx..1m: Saturated Model with Cov Matrices & Means
Mx.EM1m <-M3.1; Mx.EC1m <-X3.1; Mx.LL1m <- F3.1;

mymatrices3 <- runmx(mx.script="univSatR2.mx",mx.location="/usr/local/bin/mxt.169")
attach(mymatrices3) #matrixName groupNumber . jobNumber
#example Mx..2: Saturated Model with Raw Data
Mx.EM2 <-M3.1; Mx.EC2 <-X3.1; Mx.LL2 <- F3.1;

#OpenMx summary
cov <- rbind(cbind(EC1,EC2),cbind(EC3,EC4)); mean <- rbind(EM2,EM4); like <- rbind(cbind(LL1,LL2),cbind(LL3,LL4))
cov; mean; like

#old Mx summary
Mx.cov <- cbind(Mx.EC1,Mx.EC2); Mx.mean <- cbind(Mx.EM1m,Mx.EM2); Mx.like <- cbind(Mx.LL1,Mx.LL2)
Mx.cov; Mx.mean; Mx.like


#Compare OpenMx results to Mx results (LL: likelihood; EC: expected covariance, EM: expected means)
#1:CovPat
omxCheckCloseEnough(LL1,Mx.LL1,.001)
omxCheckCloseEnough(EC1,Mx.EC1,.001)
#2:RawPat 
omxCheckCloseEnough(LL2,Mx.LL2,.001)
omxCheckCloseEnough(EC2,Mx.EC2,.001)
omxCheckCloseEnough(EM2,Mx.EM2,.001)
#3:CovMat
omxCheckCloseEnough(LL3,Mx.LL1,.001)
omxCheckCloseEnough(EC3,Mx.EC1,.001)
#4:RawMat
omxCheckCloseEnough(LL4,Mx.LL2,.001)
omxCheckCloseEnough(EC4,Mx.EC2,.001)
omxCheckCloseEnough(EM4,Mx.EM2,.001)
