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
#dimnames=list(c("x","y"), c("x","y"))),
dat <- MASS::mvrnorm(n=500,mu=c(0,0),Sigma=matrix(c(1,0.5,0.5,1),2,2),empirical=T)
colnames(dat) <- c("x","y")

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
	firstTerm <- -0.5*(N)*sum(diag(SigmaInv%*%Der))
	if(verbose){message(paste0("firstTerm: ",firstTerm))}
	secondTerm <- 0.5*(N)*sum(diag((N-1)/N*SigmaInv%*%C%*%SigmaInv%*%Der))
	if(verbose){message(paste0("secondTerm: ",secondTerm))}
	Nu <- c(1.0,1.5) - mxGetExpected(m,"means")
	if(verbose){print("Nu:");print(Nu)}
	dNu_dtheta <- -c(0,mxEval(m1,m,T)[1])
	if(verbose){print("dNu_dtheta:");print(dNu_dtheta)}
	SigmaInvNu <- SigmaInv %*% t(Nu)
	thirdTerm <- -0.5*N*2*dNu_dtheta%*%SigmaInvNu
	if(verbose){message(paste0("thirdTerm: ",thirdTerm))}
	fourthTerm <- 0.5*N*Nu%*%SigmaInv%*%Der%*%SigmaInv%*%t(Nu)
	if(verbose){message(paste0("fourthTerm: ",fourthTerm))}
	return(-2*(firstTerm+secondTerm+thirdTerm+fourthTerm))
}

# No A paths, all params at MLE ####

plan <- mxComputeSequence(list(mxComputeOnce("fitfunction",c("fit","gradient","hessian")),mxComputeReportDeriv(),mxComputeReportExpectation()))
mxOption(NULL,"Analytic gradients","Yes"); mxOption(NULL,"Analytic RAM derivatives","Yes")
m0a <- mxModel(
	"TwoByTwo",
	type="RAM",
	plan,
	manifestVars = c("x","y"),
	mxPath(from=c("x","y"),arrows=2,values=0.998,labels=c("v1","v2"),free=T),
	mxPath(from="x",to="y",arrows=2,values=0.499,labels="c",free=T),
	mxPath(from="one",to=c("x","y"),arrows=1,labels=c("m1","m2"),values=0.0,free=T),
	mxData(observed=dat, type="raw")
)
m0a <- mxRun(m0a)

plan3 <- mxComputeSequence(list(mxComputeNumericDeriv(checkGradient=F,hessian=T),mxComputeReportDeriv(),mxComputeReportExpectation()))
mxOption(NULL,"Analytic gradients","No"); mxOption(NULL,"Analytic RAM derivatives","No")
m0n <- mxModel(
	"TwoByTwo",
	type="RAM",
	plan3,
	manifestVars = c("x","y"),
	mxPath(from=c("x","y"),arrows=2,values=0.998,labels=c("v1","v2"),free=T),
	mxPath(from="x",to="y",arrows=2,values=0.499,labels="c",free=T),
	mxPath(from="one",to=c("x","y"),arrows=1,labels=c("m1","m2"),values=0.0,free=T),
	mxData(observed=dat, type="raw")
)
m0n <- mxRun(m0n)
m0a$output$gradient
omxCheckCloseEnough(m0a$output$gradient,rep(0,5),5e-11)
m0n$output$gradient
omxCheckCloseEnough(m0n$output$gradient,rep(0,5),5e-8)
omxCheckCloseEnough(vech(m0a$output$hessian),vech(m0n$output$hessian),0.1)

# No A paths, no params at MLE ####

plan <- mxComputeSequence(list(mxComputeOnce("fitfunction",c("fit","gradient","hessian")),mxComputeReportDeriv(),mxComputeReportExpectation()))
mxOption(NULL,"Analytic gradients","Yes"); mxOption(NULL,"Analytic RAM derivatives","Yes")
m0a <- mxModel(
	"TwoByTwo",
	type="RAM",
	plan,
	manifestVars = c("x","y"),
	mxPath(from=c("x","y"),arrows=2,values=0.11,labels=c("v1","v2"),free=T),
	mxPath(from="x",to="y",arrows=2,values=0.1,labels="c",free=T),
	mxPath(from="one",to=c("x","y"),arrows=1,labels=c("m1","m2"),values=0.2,free=T),
	mxData(observed=dat, type="raw")
)
m0a <- mxRun(m0a)

plan3 <- mxComputeSequence(list(mxComputeNumericDeriv(checkGradient=F,hessian=T),mxComputeReportDeriv(),mxComputeReportExpectation()))
mxOption(NULL,"Analytic gradients","No"); mxOption(NULL,"Analytic RAM derivatives","No")
m0n <- mxModel(
	"TwoByTwo",
	type="RAM",
	plan3,
	manifestVars = c("x","y"),
	mxPath(from=c("x","y"),arrows=2,values=0.11,labels=c("v1","v2"),free=T),
	mxPath(from="x",to="y",arrows=2,values=0.1,labels="c",free=T),
	mxPath(from="one",to=c("x","y"),arrows=1,labels=c("m1","m2"),values=0.2,free=T),
	mxData(observed=dat, type="raw")
)
m0n <- mxRun(m0n)
m0a$output$gradient
m0n$output$gradient
omxCheckCloseEnough(m0a$output$gradient,m0n$output$gradient,1e-4)

# Only S paths free ####

plan <- mxComputeSequence(list(mxComputeOnce("fitfunction",c("fit","gradient","hessian")),mxComputeReportDeriv(),mxComputeReportExpectation()))
mxOption(NULL,"Analytic gradients","Yes"); mxOption(NULL,"Analytic RAM derivatives","Yes")
m1a <- mxModel(
	"TwoByTwo",
	type="RAM",
	plan,
	manifestVars = c("x","y"),
	mxPath(from=c("x","y"),arrows=2, values=0.1,labels=c("v1","v2"),free=T),
	mxPath(from="x",to="y",arrows=1,values=0.5,labels="b",free=F),
	mxPath(from="one",to=c("x","y"),arrows=1,labels=c("m1","m2"),values=1.0,free=F),
	mxData(observed=dat, type="raw")
)
m1a <- mxRun(m1a)
m1a$output$gradient
#gb(m1a,T)

plan3 <- mxComputeSequence(list(mxComputeNumericDeriv(checkGradient=F,hessian=T),mxComputeReportDeriv(),mxComputeReportExpectation()))
mxOption(NULL,"Analytic gradients","No"); mxOption(NULL,"Analytic RAM derivatives","No")
m1n <- mxModel(
	"TwoByTwo",
	type="RAM",
	plan3,
	manifestVars = c("x","y"),
	mxPath(from=c("x","y"),arrows=2, values=0.1,labels=c("v1","v2"),free=T),
	mxPath(from="x",to="y",arrows=1,values=0.5,labels="b",free=F),
	mxPath(from="one",to=c("x","y"),arrows=1,labels=c("m1","m2"),values=1.0,free=F),
	mxData(observed=dat, type="raw")
)
m1n <- mxRun(m1n)
m1a$output$gradient
m1n$output$gradient

omxCheckCloseEnough(m1a$output$gradient-m1n$output$gradient, c(0,0), 5e-7)


# Only A paths free ####

mxOption(NULL,"Analytic gradients","Yes"); mxOption(NULL,"Analytic RAM derivatives","Yes")
m2a <- mxModel(
	"TwoByTwo",
	type="RAM",
	plan,
	manifestVars = c("x","y"),
	mxPath(from=c("x","y"),arrows=2, values=c(0.998,0.7485),labels=c("v1","v2"),free=F),
	mxPath(from="x",to="y",arrows=1,values=0.1,labels="b",free=T),
	mxPath(from="one",to=c("x","y"),arrows=1,labels=c("m1","m2"),values=1.0,free=F),
	mxData(observed=dat, type="raw")
)
m2a <- mxRun(m2a)
gb(m2a,T)
m2a$output$gradient

mxOption(NULL,"Analytic gradients","No"); mxOption(NULL,"Analytic RAM derivatives","No")
m2n <- mxModel(
	"TwoByTwo",
	type="RAM",
	plan3,
	manifestVars = c("x","y"),
	mxPath(from=c("x","y"),arrows=2, values=c(0.998,0.7485),labels=c("v1","v2"),free=F),
	mxPath(from="x",to="y",arrows=1,values=0.1,labels="b",free=T),
	mxPath(from="one",to=c("x","y"),arrows=1,labels=c("m1","m2"),values=1.0,free=F),
	mxData(observed=dat, type="raw")
)
m2n <- mxRun(m2n)
m2n$output$gradient

omxCheckCloseEnough(m2a$output$gradient-m2n$output$gradient, 0, 5e-7)


# Only M paths free ####

mxOption(NULL,"Analytic gradients","Yes"); mxOption(NULL,"Analytic RAM derivatives","Yes")
m3a <- mxModel(
	"TwoByTwo",
	type="RAM",
	plan,
	manifestVars = c("x","y"),
	mxPath(from=c("x","y"),arrows=2, values=c(0.998,0.7485),labels=c("v1","v2"),free=F),
	mxPath(from="x",to="y",arrows=1,values=0.5,labels="b",free=F),
	mxPath(from="one",to=c("x","y"),arrows=1,labels=c("m1","m2"),values=0.1,free=T),
	mxData(observed=dat, type="raw")
)
m3a <- mxRun(m3a)
m3a$output$gradient

mxOption(NULL,"Analytic gradients","No"); mxOption(NULL,"Analytic RAM derivatives","No")
m3n <- mxModel(
	"TwoByTwo",
	type="RAM",
	plan3,
	manifestVars = c("x","y"),
	mxPath(from=c("x","y"),arrows=2, values=c(0.998,0.7485),labels=c("v1","v2"),free=F),
	mxPath(from="x",to="y",arrows=1,values=0.5,labels="b",free=F),
	mxPath(from="one",to=c("x","y"),arrows=1,labels=c("m1","m2"),values=0.1,free=T),
	mxData(observed=dat, type="raw")
)
m3n <- mxRun(m3n)
m3n$output$gradient

omxCheckCloseEnough(m3a$output$gradient-m3n$output$gradient, c(0,0), 5e-7)


# A & M paths free ####

mxOption(NULL,"Analytic gradients","Yes"); mxOption(NULL,"Analytic RAM derivatives","Yes")
m4a <- mxModel(
	"TwoByTwo",
	type="RAM",
	plan,
	manifestVars = c("x","y"),
	mxPath(from=c("x","y"),arrows=2, values=c(0.998,0.7485),labels=c("v1","v2"),free=F),
	mxPath(from="x",to="y",arrows=1,values=0.1,labels="b",free=T),
	mxPath(from="one",to=c("x","y"),arrows=1,labels=c("m1","m2"),values=0.2,free=T),
	mxData(observed=dat, type="raw")
)
m4a <- mxRun(m4a)
m4a$output$gradient

mxOption(NULL,"Analytic gradients","No"); mxOption(NULL,"Analytic RAM derivatives","No")
m4n <- mxModel(
	"TwoByTwo",
	type="RAM",
	plan3,
	manifestVars = c("x","y"),
	mxPath(from=c("x","y"),arrows=2, values=c(0.998,0.7485),labels=c("v1","v2"),free=F),
	mxPath(from="x",to="y",arrows=1,values=0.1,labels="b",free=T),
	mxPath(from="one",to=c("x","y"),arrows=1,labels=c("m1","m2"),values=0.2,free=T),
	mxData(observed=dat, type="raw")
)
m4n <- mxRun(m4n)
m4n$output$gradient

omxCheckCloseEnough(m4a$output$gradient-m4n$output$gradient, c(0,0,0), 5e-7)


## A paths & m2 free ####

mxOption(NULL,"Analytic gradients","Yes"); mxOption(NULL,"Analytic RAM derivatives","Yes")
m5a <- mxModel(
	"TwoByTwo",
	type="RAM",
	plan,
	manifestVars = c("x","y"),
	mxPath(from=c("x","y"),arrows=2, values=c(0.998,0.7485),labels=c("v1","v2"),free=F),
	mxPath(from="x",to="y",arrows=1,values=0.1,labels="b",free=T),
	mxPath(from="one",to=c("x","y"),arrows=1,labels=c("m1","m2"),values=c(1.0,0.1),free=c(F,T)),
	mxData(observed=dat, type="raw")
)
m5a <- mxRun(m5a)
m5a$output$gradient

mxOption(NULL,"Analytic gradients","No"); mxOption(NULL,"Analytic RAM derivatives","No")
m5n <- mxModel(
	"TwoByTwo",
	type="RAM",
	plan3,
	manifestVars = c("x","y"),
	mxPath(from=c("x","y"),arrows=2, values=c(0.998,0.7485),labels=c("v1","v2"),free=F),
	mxPath(from="x",to="y",arrows=1,values=0.1,labels="b",free=T),
	mxPath(from="one",to=c("x","y"),arrows=1,labels=c("m1","m2"),values=c(1.0,0.1),free=c(F,T)),
	mxData(observed=dat, type="raw")
)
m5n <- mxRun(m5n)
m5n$output$gradient

omxCheckCloseEnough(m5a$output$gradient-m5n$output$gradient, c(0,0), 5e-8)


## A paths & m1 free ####

mxOption(NULL,"Analytic gradients","Yes"); mxOption(NULL,"Analytic RAM derivatives","Yes")
m6a <- mxModel(
	"TwoByTwo",
	type="RAM",
	plan,
	manifestVars = c("x","y"),
	mxPath(from=c("x","y"),arrows=2, values=c(0.998,0.7485),labels=c("v1","v2"),free=F),
	mxPath(from="x",to="y",arrows=1,values=0.1,labels="b",free=T),
	mxPath(from="one",to=c("x","y"),arrows=1,labels=c("m1","m2"),values=c(0.2,1.0),free=c(T,F)),
	mxData(observed=dat, type="raw")
)
m6a <- mxRun(m6a)
m6a$output$gradient

mxOption(NULL,"Analytic gradients","No"); mxOption(NULL,"Analytic RAM derivatives","No")
m6n <- mxModel(
	"TwoByTwo",
	type="RAM",
	plan3,
	manifestVars = c("x","y"),
	mxPath(from=c("x","y"),arrows=2, values=c(0.998,0.7485),labels=c("v1","v2"),free=F),
	mxPath(from="x",to="y",arrows=1,values=0.1,labels="b",free=T),
	mxPath(from="one",to=c("x","y"),arrows=1,labels=c("m1","m2"),values=c(0.2,1.0),free=c(T,F)),
	mxData(observed=dat, type="raw")
)
m6n <- mxRun(m6n)
m6n$output$gradient

omxCheckCloseEnough(m6a$output$gradient-m6n$output$gradient, c(0,0), 5e-8)


# All paths free ####

mxOption(NULL,"Analytic gradients","Yes"); mxOption(NULL,"Analytic RAM derivatives","Yes")
m7a <- mxModel(
	"TwoByTwo",
	type="RAM",
	plan,
	manifestVars = c("x","y"),
	mxPath(from=c("x","y"),arrows=2, values=0.5,labels=c("v1","v2"),free=T),
	mxPath(from="x",to="y",arrows=1,values=0.1,labels="b",free=T),
	mxPath(from="one",to=c("x","y"),arrows=1,labels=c("m1","m2"),values=0.2,free=T),
	mxData(observed=dat, type="raw")
)
m7a <- mxRun(m7a)
m7a$output$gradient
m7a$output$hessian

mxOption(NULL,"Analytic gradients","No"); mxOption(NULL,"Analytic RAM derivatives","No")
m7n <- mxModel(
	"TwoByTwo",
	type="RAM",
	plan3,
	manifestVars = c("x","y"),
	mxPath(from=c("x","y"),arrows=2, values=0.5,labels=c("v1","v2"),free=T),
	mxPath(from="x",to="y",arrows=1,values=0.1,labels="b",free=T),
	mxPath(from="one",to=c("x","y"),arrows=1,labels=c("m1","m2"),values=0.2,free=T),
	mxData(observed=dat, type="raw")
)
m7n <- mxRun(m7n)
m7n$output$gradient
m7n$output$hessian

omxCheckCloseEnough(m7a$output$gradient-m7n$output$gradient, rep(0,5), 5e-7)

mxOption(reset=TRUE)
