require(OpenMx)
require(MASS)
setwd('/Users/hermine/Applications/bin/OpenMx/trunk/demo')
set.seed(200); rs=.5; xy <- mvrnorm (1000, c(0,0), matrix(c(1,rs,rs,1),2,2))
testData <- xy; selVars <- c('X','Y'); dimnames(testData) <- list(NULL, selVars)
summary(testData); cov(testData)

#example 1: Saturated Model with Cov Matrices and Paths Input
multSatModel1 <- mxModel("multSat1",
	manifestVars= selVars ,
	mxPath(from="X", to="Y", arrows=2, free=T, values=.2, lbound=.01, labels="covXY"),
	mxPath(from=c("X", "Y"), arrows=2, free=T, values=1, lbound=.01, labels=c("varX","varY")),
	mxData(cov(testData), type="cov", numObs=1000), type="RAM")
#	mxData(cov(testData), type="cov", numObs=1000, means=colMeans(testData)), type="RAM")
multSatFit1 <- mxRun(multSatModel1)
EC1 <- multSatFit1[['S']]@values; LL1 <- mxEvaluate(objective,multSatFit1);

#examples 2: Saturated Model with Raw Data and Paths Input
multSatModel2 <- mxModel("multSat2",
	manifestVars= selVars,
	mxPath(from="X", to="Y", arrows=2, free=T, values=.2, lbound=.01, labels="covXY"),
	mxPath(from=c("X", "Y"), arrows=2, free=T, values=1, lbound=.01, labels=c("varX","varY")),
	mxData(testData, type="raw"), type="RAM")
multSatFit2 <- mxRun(multSatModel2)
EM2 <- multSatFit2[['M']]@values; EC2 <- multSatFit2[['S']]@values; LL2 <- mxEvaluate(objective,multSatFit2);

#example 3: Saturated Model with Cov Matrices and Matrices Input
multSatModel3 <- mxModel("multSat3",
 	mxMatrix("Symm", nrow=2, ncol=2, free=T, values=c(1,.5,1), name="expCov"),
 	mxData(cov(testData), type="cov", numObs=1000),
 	mxMLObjective("expCov"))
multSatFit3 <- mxRun(multSatModel3)
EC3 <- multSatFit3[['expCov']]@values; LL3 <- mxEvaluate(objective,multSatFit3);

#examples 4: Saturated Model with Raw Data and Matrices Input
multSatModel4 <- mxModel("multSat4",
 	mxMatrix("Symm", nrow=2, ncol=2, free=T, values=c(1,.5,1), name="expCov", dimnames=list(selVars, selVars)),
 	mxMatrix("Full", nrow=1, ncol=2, free=T, values=c(0,0), name="expMean", dimnames=list(NULL, selVars)),
 	mxData(testData, type="raw"),
 	mxFIMLObjective("expCov", "expMean"))
multSatFit4 <- mxRun(multSatModel4)
EM4 <- multSatFit4[['expMean']]@values; EC4 <- multSatFit4[['expCov']]@values; LL4 <- mxEvaluate(objective,multSatFit4);

cov <- rbind(cbind(EC1,EC2),cbind(EC3,EC4)); mean <- rbind(EM2,EM4); like <- rbind(cbind(LL1,LL2),cbind(LL3,LL4))
cov; mean; like

setwd("/Users/hermine/Applications/bin/OpenMx/trunk/demo/MxR")
cat("*\n",file='star')
write.table(cov(testData),file="cov",row.names=F,quote=F,col.names=F); system("cat star cov > testData.cov")
write.table(colMeans(testData),file="mea",row.names=F,quote=F,col.names=F); system("cat star mea > testData.mea")
testDataDF<-as.data.frame(testData)
write.table(testDataDF,file="testData.rec",row.names=F,na=".",quote=F,col.names=F)

source("runmx.R")
#Mx file includes 3 Saturated Model jobs: 1: Cov Matrices, 2: Cov Matrices + Means, 3: Raw Data 
mymatrices <- runmx(mx.script="bivSatR.mx",mx.location="/usr/local/bin/mxt.169")
attach(mymatrices) #matrixName groupNumber . jobNumber
#example Mx..1: Saturated Model with Cov Matrices
Mx.EC1 <-X3.1; Mx.LL1 <- F3.1;
#example Mx..1m: Saturated Model with Cov Matrices & Means
Mx.EM1m <-M3.2; Mx.EC1m <-X3.2; Mx.LL1 <- F3.2;
#example Mx..2: Saturated Model with Raw Data
Mx.EM2 <-M3.3; Mx.EC2 <-X3.3; Mx.LL1 <- F3.3;

Mx.cov <- cbind(Mx.EC1,Mx.EC2); Mx.mean <- cbind(Mx.EM1m,Mx.EM2); Mx.like <- cbind(Mx.LL1,Mx.LL2)
Mx.cov; Mx.mean; Mx.like

#Compare OpenMx results to Mx results (LL: likelihood; EC: expected covariance, EM: expected means)
#1:CovPat
omxCheckCloseEnough(LL1,Mx.LL1,.001); omxCheckCloseEnough(EC1,Mx.EC1,.001);
#2:RawPat 
omxCheckCloseEnough(LL2,Mx.LL2,.001); omxCheckCloseEnough(EC2,Mx.EC2,.001); omxCheckCloseEnough(EM2,Mx.EM2,.001)
#3:CovMat
omxCheckCloseEnough(LL3,Mx.LL1,.001) ;omxCheckCloseEnough(EC3,Mx.EC1,.001);
#4:RawMat
omxCheckCloseEnough(LL4,Mx.LL2,.001); omxCheckCloseEnough(EC4,Mx.EC2,.001); omxCheckCloseEnough(EM4,Mx.EM2,.001)
