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
x <- rnorm(500)
x <- (x-mean(x))/sd(x)

gv <- function(m,verbose=TRUE,N=500){
	Sigma <- mxGetExpected(m,"covariance")
	C <- matrix(1,1,1)
	SigmaInv <- solve(Sigma)
	firstTerm <- -0.5*(N)*sum(diag(SigmaInv%*%matrix(1,1,1)))
	if(verbose){message(paste0("firstTerm: ",firstTerm))}
	#secondTerm <- sum(diag((N-1)/N*SigmaInv%*%C%*%SigmaInv%*%matrix(1,1,1)))
	secondTerm <- 0.5*(N)*sum(diag((N-1)/N*SigmaInv%*%C%*%SigmaInv%*%matrix(1,1,1)))
	if(verbose){message(paste0("secondTerm: ",secondTerm))}
	obmean <- mxGetExpected(m,"mean")
	Nu <- 0 - obmean
	dNu_dv <- 0
	thirdTerm <- -0.5*N*2*dNu_dv*SigmaInv*Nu
	if(verbose){message(paste0("thirdTerm: ",thirdTerm))}
	fourthTerm <- 0.5*N*Nu*SigmaInv%*%matrix(1,1,1)%*%SigmaInv*Nu
	#fourthTerm <- -0.5*(N-1)*Nu*SigmaInv%*%matrix(1,1,1)%*%SigmaInv*Nu
	if(verbose){message(paste0("fourthTerm: ",fourthTerm))}
	return(-2*(firstTerm+secondTerm+thirdTerm+fourthTerm))
}

# Check to see if gradient is zero at MLE ####

plan <- mxComputeSequence(list(mxComputeOnce("fitfunction",c("fit","gradient")),mxComputeReportDeriv(),mxComputeReportExpectation()))
mxOption(NULL,"Analytic gradients","Yes")
m1a <- mxModel(
	"Simple",
	type="RAM",
	plan,
	manifestVars = "x",
	mxPath(from="x", to="x",arrows=2, values=0.998,labels="v",free=T),
	mxPath(from="one", to="x", arrows=1, values=0.0, labels="a",free=T),
	mxData(matrix(x,dimnames=list(NULL,"x")), type="raw")
)
m1a <- mxRun(m1a)

plan3 <- mxComputeSequence(list(mxComputeNumericDeriv(checkGradient=F,hessian=F),mxComputeReportDeriv(),mxComputeReportExpectation()))
mxOption(NULL,"Analytic gradients","No")
m1n <- mxModel(
	"Simple",
	type="RAM",
	plan3,
	manifestVars = "x",
	mxPath(from="x", to="x",arrows=2, values=0.998,labels="v",free=T),
	mxPath(from="one", to="x", arrows=1, values=0.0, labels="a",free=T),
	mxData(matrix(x,dimnames=list(NULL,"x")), type="raw")
)
m1n <- mxRun(m1n)
#Make sure analytic gradient w/r/t variance is zero at the MLE in backend & frontend:
omxCheckCloseEnough(gv(m1a)[1],0,1e-12)
omxCheckCloseEnough(0,m1a$output$gradient[1],1e-12)
#More importantly, make sure analytic & numeric gradients are both zero at the MLE:
omxCheckCloseEnough(m1n$output$gradient,c(0,0),2e-8)
omxCheckCloseEnough(m1a$output$gradient,c(0,0),1e-8)

# Check to see if analytic & numeric gradients match when not at MLE ####
plan <- mxComputeSequence(list(mxComputeOnce("fitfunction",c("fit","gradient")),mxComputeReportDeriv(),mxComputeReportExpectation()))
mxOption(NULL,"Analytic gradients","Yes")
m1a <- mxModel(
	"Simple",
	type="RAM",
	plan,
	manifestVars = "x",
	mxPath(from="x", to="x",arrows=2, values=0.5,labels="v",free=T),
	mxPath(from="one", to="x", arrows=1, values=0.2, labels="a",free=T),
	mxData(matrix(x,dimnames=list(NULL,"x")), type="raw")
)
m1a <- mxRun(m1a)

plan3 <- mxComputeSequence(list(mxComputeNumericDeriv(checkGradient=F,hessian=F),mxComputeReportDeriv(),mxComputeReportExpectation()))
mxOption(NULL,"Analytic gradients","No")
m1n <- mxModel(
	"Simple",
	type="RAM",
	plan3,
	manifestVars = "x",
	mxPath(from="x", to="x",arrows=2, values=0.5,labels="v",free=T),
	mxPath(from="one", to="x", arrows=1, values=0.2, labels="a",free=T),
	mxData(matrix(x,dimnames=list(NULL,"x")), type="raw")
)
m1n <- mxRun(m1n)
#Make sure analytic gradient w/r/t variance matches in backend & frontend:
omxCheckCloseEnough(gv(m1a)[1],m1a$output$gradient[1],5e-12)
#More importantly, make sure analytic & numeric gradients match:
omxCheckCloseEnough(m1a$output$gradient,m1n$output$gradient,2e-8)
