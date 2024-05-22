#
#   Copyright 2007-2023 by the individuals mentioned in the source code history
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

require(OpenMx)
#This script does not actually invoke the optimizer, so there's no need to test it with all 3:
if(mxOption(NULL,"Default optimizer")!="SLSQP"){stop("SKIP")} 
require(numDeriv)
g <- function(m,verbose=FALSE,N=500){
	Sigma <- mxGetExpected(m,"covariance")
	C <- matrix(1,1,1)
	SigmaInv <- solve(Sigma)
	CinvDer_trace <- sum(diag(SigmaInv%*%matrix(1,1,1)))
	if(verbose){message(paste0("CinvDer_trace: ",CinvDer_trace))}
	secondTerm <- sum(diag((N-1)/N*SigmaInv%*%C%*%SigmaInv%*%matrix(1,1,1)))
	if(verbose){message(paste0("secondTerm: ",secondTerm))}
	return((N)*(CinvDer_trace - secondTerm))
}
mxf <- function(x,md){
	md@compute <- NULL
	md <- omxSetParameters(model=md,labels="v",values=x)
	mm <- mxRun(md,useOptimizer=F)
	return(mm$output$fit)
}
truevals <- 
	c(0.1,0.17349,0.21204,0.21030,0.23219,0.35802,0.55414,0.77777,0.8,0.998,1.0,1.1,1.17349,1.21030,1.21204,1.23219,1.35802,1.55414,1.77777,1.8,1.998,2.0)
results <- 
	data.frame(
		truevals=truevals,m1aBack=rep(NA_real_,22),m1aFront=rep(NA_real_,22),m2aBack=rep(NA_real_,22),m2aFront=rep(NA_real_,22),
		m3aBack=rep(NA_real_,22),m3aFront=rep(NA_real_,22),m3aFront2=rep(NA_real_,22),m3aFront3=rep(NA_real_,22)
	)
wallTimes <- data.frame(
	truevals=truevals, m1a=rep(NA,22), m2a=rep(NA,22)
)

for(i in 1:22){
	
	
	plan <- mxComputeSequence(list(mxComputeOnce("fitfunction",c("fit","gradient")),mxComputeReportDeriv(),mxComputeReportExpectation()))
	mxOption(NULL,"Analytic gradients","Yes")
	m1 <- mxModel(
		"Simple",
		type="RAM",
		plan,
		manifestVars = "x",
		mxPath(from="x", to="x",arrows=2, values=truevals[i],labels="v"),
		mxData(matrix(data=1,nrow=1,ncol=1,dimnames=list("x","x")), type="cov", numObs=500)
	)
	m1a <- mxRun(m1)
	results$m1aBack[i] <- m1a$output$gradient[1]
	results$m1aFront[i] <- g(m1a)[1]
	coef(m1a)
	wallTimes$m1a[i] <- summary(m1a)$wallTime
	
	
	mxOption(NULL,"Analytic gradients","No")
	m2a <- mxRun(m1)
	results$m2aBack[i] <- m2a$output$gradient[1]
	results$m2aFront[i] <- g(m2a)[1]
	coef(m2a)
	wallTimes$m2a[i] <- summary(m2a)$wallTime
	
	plan3 <- mxComputeSequence(list(mxComputeNumericDeriv(checkGradient=F,hessian=F),mxComputeReportDeriv(),mxComputeReportExpectation()))
	m3 <- mxModel(
		"Simple",
		type="RAM",
		plan3,
		manifestVars = "x",
		mxPath(from="x", to="x",arrows=2, values=truevals[i],labels="v"),
		mxData(matrix(data=1,nrow=1,ncol=1,dimnames=list("x","x")), type="cov", numObs=500)
	)
	m3a <- mxRun(m3)
	results$m3aBack[i] <- m3a$output$gradient[1]
	results$m3aFront[i] <- g(m3a)[1]
	results$m3aFront3[i] <- numDeriv::grad(mxf,truevals[i],md=m3a)
	coef(m3a)
	
	print(i)
}

# Frontend and backend analytic results should agree:
omxCheckCloseEnough(cor(results$m1aBack,results$m1aFront),1.0,1e-4)
omxCheckCloseEnough(results$m1aBack-results$m1aFront,rep(0.0,22),1e-10)

# Gradient should be zero at MLE:
omxCheckCloseEnough(results$m1aBack[10],0.0,1e-10)
omxCheckCloseEnough(results$m1aFront[10],0.0,1e-10)

#OpenMx's and numDeriv's numeric results should agree:
omxCheckCloseEnough(cor(results$m3aBack, results$m3aFront3),1.0,1e-4)
omxCheckCloseEnough(results$m3aBack-results$m3aFront3,rep(0.0,22),1e-6)

#Gradient should be zero at MLE:
omxCheckCloseEnough(results$m3aBack[10],0.0,1e-8)
omxCheckCloseEnough(results$m3aFront[10],0.0,1e-8)

#Analytic and numeric results should agree:
omxCheckCloseEnough(results$m1aBack-results$m3aBack,rep(0.0,22),1e-7)
omxCheckCloseEnough((results$m3aBack/results$m1aBack)[-10],rep(1.0,21),1e-8)

#Analytic gradient should be faster:
if(0){
	omxCheckTrue(all(wallTimes$m1a < wallTimes$m2a)) #<--FALSE
}
wallTimes$m1a - wallTimes$m2a
