#
#   Copyright 2007-2018 The OpenMx Project
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

library(OpenMx)
#mxOption(NULL,"Default optimizer","SLSQP")

#CSOLNP still fails in a platform-specific manner:
if(mxOption(NULL,"Default optimizer")=="CSOLNP" && .Platform$OS.type=="windows" && .Platform$r_arch=="i386"){stop("SKIP")}

library(mvtnorm)

set.seed(170209)
Rtrue <- matrix(c(1,0.5,0.5,0.5,1,0.5,0.5,0.5,1),nrow=3,byrow=T)
latdat <- rmvnorm(n=5000,mean=rep(0,3),sigma=Rtrue)
latdat[,1] <- ifelse(latdat[,1]>1.64,1,0)
latdat[,2] <- ifelse(latdat[,2]>1.64,1,0)
latdat[,3] <- ifelse(latdat[,3]>1.64,1,0)
table(apply(latdat,1,function(x){paste(x[1],x[2],x[3],sep="")}))

obsdat <- data.frame(
	y1=mxFactor(latdat[,1],levels=c(0,1)),
	y2=mxFactor(latdat[,2],levels=c(0,1)),
	y3=mxFactor(latdat[,3],levels=c(0,1))
)

mxdat <- mxData(observed=obsdat,type="raw")

safeT <- mxConstraint(-1%x%Sigma<Zilch,name="safety")
safeT2 <- mxConstraint(-1%x%Sigma<Zilch,name="safety",jac="ineqJac")
sgn <- mxMatrix(type="Full",nrow=1,ncol=1,free=F, #This is weird but necessary.
								values=ifelse(mxOption(NULL,"Default optimizer")=="NPSOL",-1,1),
								#values=1,
								name="Sgn")

#mxOption(NULL,"Standard Errors","No")
#mxOption(NULL,"Calculate Hessian","No")
#plan <- omxDefaultComputePlan(modelName="mod1")
#plan$steps$GD$verbose <- 5L

m1 <- mxModel(
	"mod1",
	mxdat,
	#plan,
	mxMatrix(type="Lower",nrow=3,free=T,values=c(1,1e-7,1e-7,1,1e-7,1),labels=c("s11","s21","s31","s22","s32","s33"),name="S"),
	mxMatrix(type="Zero",nrow=1,ncol=3,name="Mu"),
	mxMatrix(type="Unit",nrow=3,ncol=1,name="ONE"),
	mxMatrix(type="Zero",nrow=3,ncol=3,name="Zilch"),
	mxMatrix(type="Full",nrow=1,ncol=3,free=T,values=1.64,labels=c("tau1","tau2","tau3"),name="Tau"),
	mxAlgebra(S%*%t(S),name="Sigma"),
	mxFitFunctionML(),
	mxExpectationNormal(covariance="Sigma",means="Mu",dimnames=c("y1","y2","y3"),thresholds="Tau",threshnames=c("y1","y2","y3")),
	safeT,
	mxConstraint(diag2vec(Sigma)==ONE,name="identifying")
)
m2 <- mxRun(m1)
summary(m2)
mxEval(Sigma,m2,T)
omxCheckCloseEnough(mxEval(Tau,m2,T)[1,],c(1.64,1.64,1.64),0.1)
omxCheckCloseEnough(mxEval(Sigma,m2,T)[c(2,3,6)],c(0.5,0.5,0.5),0.05)
omxCheckCloseEnough(diag(mxEval(Sigma,m2,T)),c(1,1,1),as.numeric(mxOption(NULL,"Feasibility tolerance")))



eqjsub <- mxAlgebra(rbind(
	#s11		s21		s31		s22		s32		s33
	cbind(2*S[1,1], 0,  0,    0,    0,    0),
	cbind(0, 2*S[2,1], 0, 2*S[2,2], 0, 0),
	cbind(0, 0, 2*S[3,1], 0, 2*S[3,2], 2*S[3,3])
), name="eqJacSub")
taueqjac <- mxMatrix(type="Zero",nrow=3,ncol=3,name="tauEqJac")
tauineqjac <- mxMatrix(type="Zero",nrow=9,ncol=3,name="tauIneqJac")
ineqjsub <- mxAlgebra(rbind(
	#s11		s21		s31		s22		s32		s33
	cbind(2*S[1,1],0,0,0,0,0),
	cbind(S[2,1],S[1,1],0,0,0,0),
	cbind(S[3,1],0,S[1,1],0,0,0),
	cbind(S[2,1],S[1,1],0,0,0,0),
	cbind(0,2*S[2,1],0,2*S[2,2],0,0),
	cbind(0,S[3,1],S[2,1],S[3,2],S[2,2],0),
	cbind(S[3,1],0,S[1,1],0,0,0),
	cbind(0,S[3,1],S[2,1],S[3,2],S[2,2],0),
	cbind(0,0,2*S[3,1],0,2*S[3,2],2*S[3,3])
), name="ineqJacSub")
eqjac <- mxAlgebra(cbind(eqJacSub,tauEqJac),name="eqJac",
									 dimnames=list(NULL,c("s11","s21","s31","s22","s32","s33","tau1","tau2","tau3")))
ineqjac <- mxAlgebra(Sgn%x%cbind(ineqJacSub,tauIneqJac),name="ineqJac",
										 dimnames=list(NULL,c("s11","s21","s31","s22","s32","s33","tau1","tau2","tau3")))

# mxOption(NULL,"Standard Errors","No")
# mxOption(NULL,"Calculate Hessian","No")
# plan <- omxDefaultComputePlan(modelName="mod3")
# plan$steps$GD$verbose <- 5L

m3 <- mxModel(
	"mod3",
	mxdat,
	#plan,
	mxMatrix(type="Lower",nrow=3,free=T,values=c(1,1e-7,1e-7,1,1e-7,1),labels=c("s11","s21","s31","s22","s32","s33"),name="S"),
	mxMatrix(type="Zero",nrow=1,ncol=3,name="Mu"),
	mxMatrix(type="Unit",nrow=3,ncol=1,name="ONE"),
	mxMatrix(type="Zero",nrow=3,ncol=3,name="Zilch"),
	mxMatrix(type="Full",nrow=1,ncol=3,free=T,values=1.64,labels=c("tau1","tau2","tau3"),name="Tau"),
	mxAlgebra(S%*%t(S),name="Sigma"),
	mxFitFunctionML(),
	mxExpectationNormal(covariance="Sigma",means="Mu",dimnames=c("y1","y2","y3"),thresholds="Tau",threshnames=c("y1","y2","y3")),
	safeT,#safeT2,
	mxConstraint(diag2vec(Sigma)==ONE,name="identifying",jac="eqJac"),
	eqjsub, taueqjac, eqjac#, 
	#tauineqjac, ineqjsub, ineqjac, sgn
)
m4 <- mxRun(m3)
summary(m4)
mxEval(Sigma,m4,T)
if(mxOption(NULL,"Default optimizer") %in% c("NPSOL","SLSQP")){
	omxCheckTrue(m2$output$evaluations > m4$output$evaluations)
}
omxCheckCloseEnough(mxEval(Tau,m4,T)[1,],c(1.64,1.64,1.64),0.1)
omxCheckCloseEnough(mxEval(Sigma,m4,T)[c(2,3,6)],c(0.5,0.5,0.5),0.05)
omxCheckCloseEnough(diag(mxEval(Sigma,m4,T)),c(1,1,1),as.numeric(mxOption(NULL,"Feasibility tolerance")))
omxCheckCloseEnough(mxEval(Tau,m4,T)[1,],mxEval(Tau,m2,T)[1,],5e-4)
omxCheckCloseEnough(mxEval(Sigma,m4,T),mxEval(Sigma,m2,T),5e-4)



# plan <- omxDefaultComputePlan(modelName="mod5")
# plan$steps$GD$verbose <- 5L

m5 <- mxModel(
	"mod5",
	mxdat,
	#plan,
	mxMatrix(type="Lower",nrow=3,free=T,values=c(1,1e-7,1e-7,1,1e-7,1),labels=c("s11","s21","s31","s22","s32","s33"),name="S"),
	mxMatrix(type="Zero",nrow=1,ncol=3,name="Mu"),
	mxMatrix(type="Unit",nrow=3,ncol=1,name="ONE"),
	mxMatrix(type="Zero",nrow=3,ncol=3,name="Zilch"),
	mxMatrix(type="Full",nrow=1,ncol=3,free=T,values=1.64,labels=c("tau1","tau2","tau3"),name="Tau"),
	mxAlgebra(S%*%t(S),name="Sigma"),
	mxFitFunctionML(),
	mxExpectationNormal(covariance="Sigma",means="Mu",dimnames=c("y1","y2","y3"),thresholds="Tau",threshnames=c("y1","y2","y3")),
	safeT2,
	mxConstraint(diag2vec(Sigma)==ONE,name="identifying",jac="eqJac"),
	eqjsub, taueqjac, eqjac, 
	tauineqjac, ineqjsub, ineqjac, sgn
)
m6 <- mxRun(m5)
summary(m6)
mxEval(Sigma,m6,T)
#Interestingly, SLSQP doesn't gain any advantage in function evaluations by adding analytic derivatives
#for the inequality constraints, but NPSOL does:
if(mxOption(NULL,"Default optimizer") %in% c("CSOLNP","NPSOL")){
	omxCheckTrue(m4$output$evaluations > m6$output$evaluations)
}
omxCheckCloseEnough(mxEval(Tau,m6,T)[1,],c(1.64,1.64,1.64),0.1)
omxCheckCloseEnough(mxEval(Sigma,m6,T)[c(2,3,6)],c(0.5,0.5,0.5),0.05)
omxCheckCloseEnough(diag(mxEval(Sigma,m6,T)),c(1,1,1),as.numeric(mxOption(NULL,"Feasibility tolerance")))
omxCheckCloseEnough(mxEval(Tau,m6,T)[1,],mxEval(Tau,m2,T)[1,],5e-4)
omxCheckCloseEnough(mxEval(Sigma,m6,T),mxEval(Sigma,m2,T),5e-4)
