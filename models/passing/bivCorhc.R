setwd("~/Applications/bin/OpenMx/trunk/demo/ExamplesH")
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

#Fit Saturated Model with RawData and Matrices Input using Cholesky Decomposition
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
EM <- bivCorFit[['expMean']]@values
EC <- mxEvaluate(Chol %*% t(Chol),bivCorFit)
LL <- mxEvaluate(objective,bivCorFit);
#eM <- mxEvaluate('expMean',bivCorFit)


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
EMs <- bivCorFitSub[['expMean']]@values
ECs <- mxEvaluate(Chol %*% t(Chol),bivCorFitSub)
LLs <- mxEvaluate(objective,bivCorFitSub);
Chi= LLs-LL;
LRT= rbind(LL,LLs,Chi); LRT


#Mx answers hard-coded
#example Mx..1m: Saturated Model with Raw Data
Mx.EM <- matrix(c(0.03211656, -0.004883885),1,2)
Mx.EC <- matrix(c( 1.0092853, 0.4813504, 0.4813504, 0.9935390),2,2)
Mx.LL <- 5415.772

#Compare OpenMx results to Mx results (LL: likelihood; EC: expected covariance, EM: expected means)
omxCheckCloseEnough(LL,Mx.LL,.001)
omxCheckCloseEnough(EC,Mx.EC,.001)
omxCheckCloseEnough(EM,Mx.EM,.001)

