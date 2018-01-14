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
if(mxOption(NULL,"Default optimizer")!="SLSQP"){stop("SKIP")}
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

#With SLSQP:
m1 <- mxModel(
	"mod1",
	mxdat,
	mxMatrix(type="Lower",nrow=3,free=T,values=c(1,0.01,0.01,1,0.01,1),labels=c("s11","s21","s31","s22","s32","s33"),name="S"),
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


#Use Nelder-Mead, trying the different methods for equality and inequality constraints: ##################

plan <- omxDefaultComputePlan()
plan$steps <- list(
	mxComputeNelderMead(eqConstraintMthd="soft",ineqConstraintMthd="soft"),
	plan$steps$RE
)
m3 <- mxModel(
	"mod3",
	mxdat,
	plan,
	mxMatrix(type="Lower",nrow=3,free=T,values=c(1,0.01,0.01,1,0.01,1),labels=c("s11","s21","s31","s22","s32","s33"),name="S"),
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
m4 <- mxRun(m3)
summary(m4)
mxEval(Sigma,m2,T)
mxEval(Sigma,m4,T)
m2$output$estimate
m4$output$estimate
m2$output$fit
m4$output$fit


plan <- omxDefaultComputePlan()
plan$steps <- list(
	mxComputeNelderMead(eqConstraintMthd="backtrack",ineqConstraintMthd="soft"),
	plan$steps$RE
)
m5 <- mxModel(
	"mod5",
	mxdat,
	plan,
	mxMatrix(type="Lower",nrow=3,free=T,values=c(1,0.01,0.01,1,0.01,1),labels=c("s11","s21","s31","s22","s32","s33"),name="S"),
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
m6 <- mxRun(m5)
summary(m6)
mxEval(Sigma,m6,T)
m2$output$estimate
m6$output$estimate
m2$output$fit
m6$output$fit


plan <- omxDefaultComputePlan()
plan$steps <- list(
	mxComputeNelderMead(eqConstraintMthd="l1p",ineqConstraintMthd="soft"),
	plan$steps$RE
)
m7 <- mxModel(
	"mod7",
	mxdat,
	plan,
	mxMatrix(type="Lower",nrow=3,free=T,values=c(1,0.01,0.01,1,0.01,1),labels=c("s11","s21","s31","s22","s32","s33"),name="S"),
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
m8 <- mxRun(m7)
summary(m8)
mxEval(Sigma,m8,T)
m2$output$estimate
m8$output$estimate
m2$output$fit
m8$output$fit
m8$compute$steps[[1]]$output$penalizedFit


plan <- omxDefaultComputePlan()
plan$steps <- list(
	mxComputeNelderMead(eqConstraintMthd="backtrack",ineqConstraintMthd="eqMthd"),
	plan$steps$RE
)
m9 <- mxModel(
	"mod9",
	mxdat,
	plan,
	mxMatrix(type="Lower",nrow=3,free=T,values=c(1,0.01,0.01,1,0.01,1),labels=c("s11","s21","s31","s22","s32","s33"),name="S"),
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
m10 <- mxRun(m9)
summary(m10)
mxEval(Sigma,m10,T)
m2$output$estimate
m10$output$estimate
m2$output$fit
m10$output$fit


plan <- omxDefaultComputePlan()
plan$steps <- list(
	mxComputeNelderMead(eqConstraintMthd="l1p",ineqConstraintMthd="eqMthd"),
	plan$steps$RE
)
m11 <- mxModel(
	"mod11",
	mxdat,
	plan,
	mxMatrix(type="Lower",nrow=3,free=T,values=c(1,0.01,0.01,1,0.01,1),labels=c("s11","s21","s31","s22","s32","s33"),name="S"),
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
m12 <- mxRun(m11)
summary(m12)
mxEval(Sigma,m12,T)
m2$output$estimate
m12$output$estimate
m2$output$fit
m12$output$fit
m12$compute$steps[[1]]$output$penalizedFit

#Holding the other Nelder-Mead arguments at their defaults, only GDsearch gets good results, and only with a stricter feasibility tolerance:

mxOption(NULL,"Feasibility tolerance",0.001)

plan <- omxDefaultComputePlan()
plan$steps <- list(
	mxComputeNelderMead(eqConstraintMthd="GDsearch",ineqConstraintMthd="soft",verbose=0L),
	plan$steps$RE
)
m13 <- mxModel(
	"mod13",
	mxdat,
	plan,
	mxMatrix(type="Lower",nrow=3,free=T,values=c(1,0.01,0.01,1,0.01,1),labels=c("s11","s21","s31","s22","s32","s33"),name="S"),
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
# m13 <- mxOption(m13,"Always Checkpoint","Yes")
# m13 <- mxOption(m13,"Checkpoint Units","minutes")
# m13 <- mxOption(m13,"Checkpoint Count",1)
m14 <- mxRun(m13)
summary(m14)
#Evaluations of the feasibility-search object count toward the total number of function evaluations:
m14$output$evaluations
mxEval(Sigma,m14,T)
m2$output$estimate
m14$output$estimate
omxCheckCloseEnough(m2$output$estimate, m14$output$estimate, 5e-4)
m2$output$fit
m14$output$fit
omxCheckCloseEnough(m2$output$fit, m14$output$fit, 5e-4)

#Changing ineqConstraintMthd from "soft" to "eqMthd" doesn't change anything in this case:
plan <- omxDefaultComputePlan()
plan$steps <- list(
	mxComputeNelderMead(eqConstraintMthd="GDsearch",ineqConstraintMthd="eqMthd",verbose=0L),
	plan$steps$RE
)
m15 <- mxModel(
	"mod15",
	mxdat,
	plan,
	mxMatrix(type="Lower",nrow=3,free=T,values=c(1,0.01,0.01,1,0.01,1),labels=c("s11","s21","s31","s22","s32","s33"),name="S"),
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
m16 <- mxRun(m15)
summary(m16)
m16$output$evaluations
mxEval(Sigma,m16,T)
m2$output$estimate
m16$output$estimate
omxCheckCloseEnough(m2$output$estimate, m16$output$estimate, 5e-4)
m2$output$fit
m16$output$fit
omxCheckCloseEnough(m2$output$fit, m16$output$fit, 1e-4)
