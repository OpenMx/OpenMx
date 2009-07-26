require(OpenMx)
require(MASS)

#Simulate Data
set.seed(200)
rs=.5
xy <- mvrnorm (1000, c(0,0), matrix(c(1,rs,rs,1),2,2))
testData <- xy
selVars <- c('X','Y')
dimnames(testData) <- list(NULL, selVars)
summary(testData)
cov(testData)

#Fit Saturated Model with RawData and Matrices Input
bivCorModel <- mxModel(
	"bivCor",
	mxMatrix(
		type="Full", 
		nrow=2, 
		ncol=2, 
		free=c(T,T,F,T), 
		values=c(1,.2,0,1), 
		dimnames=list(selVars, selVars), 
		name="Chol"), 
	mxAlgebra(
		Chol %*% t(Chol), 
		name="expCov", 
		dimnames=list(selVars, selVars)), 
	mxMatrix(
		type="Full", 
		nrow=1, 
		ncol=2, 
		free=T, 
		values=c(0,0), 
		dimnames=list(NULL, selVars), 
		name="expMean"), 
	mxData(
		testData, 
		type="raw"), 
	mxFIMLObjective(
		"expCov", 
		"expMean"))

bivCorFit <- mxRun(bivCorModel)
EM <- mxEvaluate(expMean, bivCorFit)
EC <- mxEvaluate(expCov, bivCorFit)
LL <- mxEvaluate(objective,bivCorFit);


#Test for Covariance=Zero
bivCorModelSub <-mxModel(bivCorModel,
 	mxMatrix(
		"Full", 
		nrow=2, 
		ncol=2, 
		free=c(T,F,F,T), 
		values=c(1,0,0,1), 
		name="Chol", 
		dimnames=list(selVars, selVars))
)
bivCorFitSub <- mxRun(bivCorModelSub)
EMs <- mxEvaluate(expMean, bivCorFitSub)
ECs <- mxEvaluate(expCov, bivCorFitSub)
LLs <- mxEvaluate(objective, bivCorFitSub);
Chi= LLs-LL;
LRT= rbind(LL,LLs,Chi); LRT

original.directory <- getwd()
setwd('temp-files')

#Save Data in Mx Format
testDataDF<-as.data.frame(testData)
write.table(testDataDF,file="testData.rec",row.names=F,na=".",quote=F,col.names=F)

setwd(original.directory)

#Run Mx
mymatrices <- omxOriginalMx("mx-scripts/bivCor.mx", "temp-files")
attach(mymatrices) #matrixName groupNumber . jobNumber
Mx.EM <-M3.1; Mx.EC <-X3.1; Mx.LL <- F3.1;

#Compare OpenMx results to Mx results 
#LL: likelihood; EC: expected covariance, EM: expected means
omxCheckCloseEnough(LL,Mx.LL,.001)
omxCheckCloseEnough(EC,Mx.EC,.001)
omxCheckCloseEnough(EM,Mx.EM,.001)
