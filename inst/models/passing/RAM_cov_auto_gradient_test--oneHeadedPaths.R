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
gb <- function(m,verbose=FALSE,N=500){
	Sigma <- mxGetExpected(m,"covariance")
	A <- m$A$values
	C <- matrix(c(1,.5,.5,1),2,2,dimnames=list(c("x","y"), c("x","y")))
	SigmaInv <- solve(Sigma)
	edA <- -1*matrix(c(0,1,0,0),2,2)
	I_A <- solve(diag(2) - A)
	if(verbose){print("I_A:"); print(I_A)}
	firstPart <- -1*I_A%*%edA%*%Sigma
	if(verbose){print("firstPart:"); print(firstPart)}
	secondPart <- matrix(0,2,2)
	thirdPart <- -1*Sigma%*%t(edA)%*%t(I_A)
	if(verbose){print("thirdPart:"); print(thirdPart)}
	Der <- firstPart+secondPart+thirdPart
	CinvDer_trace <- sum(diag(SigmaInv%*%Der))
	if(verbose){message(paste0("CinvDer_trace: ",CinvDer_trace))}
	secondTerm <- sum(diag((N-1)/N*SigmaInv%*%C%*%SigmaInv%*%Der))
	if(verbose){message(paste0("secondTerm: ",secondTerm))}
	return((N)*(CinvDer_trace - secondTerm))
}

plan <- mxComputeSequence(list(mxComputeOnce("fitfunction",c("fit","gradient")),mxComputeReportDeriv(),mxComputeReportExpectation()))
mxOption(NULL,"Analytic gradients","Yes"); mxOption(NULL,"Analytic RAM derivatives","Yes")
m1 <- mxModel(
	"TwoByTwo",
	type="RAM",
	plan,
	manifestVars = c("x","y"),
	mxPath(from=c("x","y"),arrows=2, values=c(0.998,0.7485),labels=c("v1","v2"),free=F),
	mxPath(from="x",to="y",arrows=1,values=0.5,labels="b",free=T),
	mxData(matrix(c(1,.5,.5,1),2,2,dimnames=list(c("x","y"), c("x","y"))), type="cov", numObs=500)
)
m1a <- mxRun(m1)
m1a$output$gradient
gb(m1a)
omxCheckCloseEnough(m1a$output$gradient,0,1e-9)
omxCheckCloseEnough(0,gb(m1a),1e-9)
summary(m1a)
mxOption(reset=TRUE)
