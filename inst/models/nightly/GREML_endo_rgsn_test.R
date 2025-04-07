#
#   Copyright 2007-2025 by the individuals mentioned in the source code history
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


library(mvtnorm)
library(Matrix)
library(OpenMx)
if(mxOption(NULL,"Default optimizer")!="CSOLNP"){stop("SKIP")}
options(mxCondenseMatrixSlots=TRUE)
# mxOption(NULL,"Number of Threads",7)
set.seed(476) 

#Number of simulees (participants):
N <- 1200

#True parameter values, for data generation:
truevals <- c(
	mu1=10,
	va1=5,
	ve1=5,
	b=1,
	mu2=1,
	va2=0.5,
	ve2=0.5
)

# GRM <- matrix(0,N,N)
# vechs(GRM) <- 0.3
# GRM <- GRM + t(GRM)
# diag(GRM) <- 1.0

msnps <- 50000
#Generate msnps SNPs.  SNPs are in linkage equilibrium, with MAF in interval [0.05,0.5]:
snps <- matrix(NA_real_,nrow=N,ncol=msnps)
for(mi in 1:msnps){
	maf <- runif(1,min=0.05,max=0.5)
	snps[,mi] <- rbinom(n=N,size=2,prob=maf)
	snps[,mi] <- (snps[,mi]-mean(snps[,mi]))/sd(snps[,mi])
	#print(mi)
}
GRM <- snps%*%t(snps) / msnps #<--#Calculate GRM from SNPs.
ev <- eigen(GRM,symmetric=T,only.values=T) #<--Eigen-decompose the GRM.
if(!all(ev$values > .Machine$double.eps)){
	GRM <- as.matrix(nearPD(GRM)$mat)
}

y1 <- as.vector(rmvnorm(n=1,mean=rep(truevals["mu1"],N),sigma=(truevals["va1"]*GRM)+diag(truevals["ve1"],N)))
y2 <- as.vector(rmvnorm(n=1,mean=rep(truevals["mu2"],N),sigma=(truevals["va2"]*GRM)+diag(truevals["ve2"],N))) + truevals["b"]*y1
widedata <- cbind(y1=y1,y2=y2)
lmod <- lm(y2~y1)

plan <- mxComputeSequence(
	steps=list(
		mxComputeNewtonRaphson(verbose=5L),
		mxComputeOnce('fitfunction', c('gradient','hessian'),verbose=5L),
		mxComputeStandardError(),
		mxComputeHessianQuality(),
		mxComputeReportDeriv(),
		mxComputeReportExpectation()
	))


m1 <- mxModel(
	"EndoRgsn",
	plan,
	
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(widedata[,"y1"])/2,labels="ve1",lbound=0.0001,name="Ve1"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(lmod$residuals)/2,labels="ve2",lbound=0.0001,name="Ve2"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(widedata[,"y1"])/2,labels="va1",name="Va1"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(lmod$residuals)/2,labels="va2",name="Va2"),
	
	mxMatrix(type="Unit",nrow=N,ncol=1,name="Uno"),
	mxMatrix(type="Symm",nrow=N,free=F,values=GRM,name="A"),
	
	mxAlgebra(A%x%Va1 + vec2diag(Uno%x%Ve1), name="V11"),
	
	mxMatrix(type="Zero",nrow=N,ncol=N,name="V12"),
	mxAlgebra(A%x%Va2 + vec2diag(Uno%x%Ve2), name="V22"),
	mxAlgebra(rbind(
		cbind(V11,V12),
		cbind(V12,V22)
	),name="V"),
	
	mxAlgebra(rbind(
		cbind(vec2diag(Uno),V12),
		cbind(V12,V12)
	),name="dV_dve1"),
	mxAlgebra(rbind(
		cbind(V12,V12),
		cbind(V12,vec2diag(Uno))
	),name="dV_dve2"),
	mxAlgebra(rbind(
		cbind(A,V12),
		cbind(V12,V12)
	),name="dV_dva1"),
	mxAlgebra(rbind(
		cbind(V12,V12),
		cbind(V12,A)
	),name="dV_dva2"),
	
	
	mxData(observed=widedata,type="raw",sort=F),
	mxExpectationGREML(V="V",yvars=c("y1","y2"),Xvars=list(character(0),c("y1")),addOnes=T,REML=F),
	mxFitFunctionGREML(dV=c(ve1="dV_dve1",ve2="dV_dve2",va1="dV_dva1",va2="dV_dva2"))
)
m1run <- mxRun(m1)
( m1sum <- summary(m1run) )


m2 <- mxModel(
	"EndoRgsn",
	
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(widedata[,"y1"])/2,labels="ve1",lbound=0.0001,name="Ve1"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(lmod$residuals)/2,labels="ve2",lbound=0.0001,name="Ve2"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(widedata[,"y1"])/2,labels="va1",name="Va1"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(lmod$residuals)/2,labels="va2",name="Va2"),
	
	mxMatrix(
		type="Full",nrow=2*N,ncol=1,free=T,
		values=c(rep(mean(widedata[,"y1"]),N), rep(lmod$coef[1],N)),
		labels=c(rep("mu1",N), rep("mu2",N)),
		name="M0"),
	mxMatrix(
		type="Full",nrow=2*N,ncol=1,free=F,
		values=c(rep(0,N), widedata[,"y1"]),
		name="M1"),
	mxMatrix(
		type="Full",nrow=1,ncol=1,free=T,
		values=lmod$coef[2],
		labels="b",name="B"),
	
	mxMatrix(type="Unit",nrow=N,ncol=1,name="Uno"),
	mxMatrix(type="Zero",nrow=N,ncol=N,name="Zilch"),
	mxMatrix(type="Zero",nrow=N,ncol=1,name="Zip"),
	mxMatrix(type="Symm",nrow=N,free=F,values=GRM,name="A"),
	
	mxAlgebra(M0 + B*M1, name="yhat"),
	
	mxAlgebra(A%x%Va1 + vec2diag(Uno%x%Ve1), name="V11"),
	mxMatrix(type="Zero",nrow=N,ncol=N,name="V12"),
	mxAlgebra(A%x%Va2 + vec2diag(Uno%x%Ve2), name="V22"),
	mxAlgebra(rbind(
		cbind(V11,V12),
		cbind(V12,V22)
	),name="V"),
	
	mxAlgebra(rbind(
		cbind(vec2diag(Uno),V12),
		cbind(V12,V12)
	),name="dV_dve1"),
	mxAlgebra(rbind(
		cbind(V12,V12),
		cbind(V12,vec2diag(Uno))
	),name="dV_dve2"),
	mxAlgebra(rbind(
		cbind(A,V12),
		cbind(V12,V12)
	),name="dV_dva1"),
	mxAlgebra(rbind(
		cbind(V12,V12),
		cbind(V12,A)
	),name="dV_dva2"),
	mxAlgebra(rbind(
		cbind(Zilch,Zilch),
		cbind(Zilch,Zilch)
	),name="dV_dmu"),
	mxAlgebra(rbind(Uno,Zip),name="dyhat_dmu1"),
	mxAlgebra(rbind(Zip,Uno),name="dyhat_dmu2"),
	mxAlgebra(rbind(Zip,Zip),name="dyhat_dV"),

	mxData(observed=widedata,type="raw",sort=F),
	mxExpectationGREML(V="V",yvars=c("y1","y2"),REML=F,yhat="yhat"),
	mxFitFunctionGREML(
		dV=c(ve1="dV_dve1",ve2="dV_dve2",va1="dV_dva1",va2="dV_dva2",mu1="dV_dmu",mu2="dV_dmu",b="dV_dmu"),
		dyhat=c(ve1="dyhat_dV",ve2="dyhat_dV",va1="dyhat_dV",va2="dyhat_dV",mu1="dyhat_dmu1",mu2="dyhat_dmu2",b="M1")
	)
)
m2run <- mxRun(m2)
summary(m2run)

omxCheckCloseEnough(m1run$output$fit,m2run$output$fit,1e-8)
omxCheckCloseEnough(coef(m1run),coef(m2run)[1:4],5e-5)
omxCheckCloseEnough(m1sum$GREMLfixeff$coeff,coef(m2run)[5:7],5e-9)
omxCheckWithinPercentError(m1run$output$standardErrors[1:4],m2run$output$standardErrors[1:4])
omxCheckCloseEnough(m1sum$GREMLfixeff$se,m2run$output$standardErrors[5:7],5e-7)
omxCheckCloseEnough(m1run$output$gradient,m2run$output$gradient[1:4],1e-5)


m3 <- mxModel(
	"EndoRgsn",
	
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(widedata[,"y1"])/2,labels="ve1",lbound=0.0001,name="Ve1"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(lmod$residuals)/2,labels="ve2",lbound=0.0001,name="Ve2"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(widedata[,"y1"])/2,labels="va1",lbound=0,name="Va1"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(lmod$residuals)/2,labels="va2",lbound=0,name="Va2"),
	
	mxMatrix(
		type="Full",nrow=2*N,ncol=1,free=T,
		values=c(rep(mean(widedata[,"y1"]),N), rep(lmod$coef[1],N)),
		labels=c(rep("mu1",N), rep("mu2",N)),
		name="M0"),
	mxMatrix(
		type="Full",nrow=2*N,ncol=1,
		free=c(rep(F,N), rep(T,N)),
		values=c(rep(0,N), rep(mean(widedata[,"y1"]),N)),
		labels=c(rep(NA,N), rep("mu1",N)),
		name="M1"),
	mxMatrix(
		type="Full",nrow=1,ncol=1,free=T,
		values=lmod$coef[2],
		labels="b",name="B"),
	
	mxMatrix(type="Unit",nrow=N,ncol=1,name="Uno"),
	mxMatrix(type="Zero",nrow=N,ncol=N,name="Zilch"),
	mxMatrix(type="Zero",nrow=N,ncol=1,name="Zip"),
	mxMatrix(type="Symm",nrow=N,free=F,values=GRM,name="A"),
	
	mxAlgebra(M0 + B*M1, name="yhat"),
	
	mxAlgebra(A%x%Va1 + vec2diag(Uno%x%Ve1), name="V11"),
	mxAlgebra(V11%x%B, name="V12"),
	mxAlgebra(V11%x%(B^2) + A%x%Va2 + vec2diag(Uno%x%Ve2), name="V22"),
	mxAlgebra(rbind(
		cbind(V11,V12),
		cbind(V12,V22)
	),name="V"),
	
	mxAlgebra(rbind(
		cbind(vec2diag(Uno),vec2diag(Uno%x%B)),
		cbind(vec2diag(Uno%x%B),vec2diag(Uno%x%(B^2)))
	),name="dV_dve1"),
	mxAlgebra(rbind(
		cbind(Zilch,Zilch),
		cbind(Zilch,vec2diag(Uno))
	),name="dV_dve2"),
	mxAlgebra(rbind(
		cbind(A,A%x%B),
		cbind(A%x%B,A%x%(B^2))
	),name="dV_dva1"),
	mxAlgebra(rbind(
		cbind(Zilch,Zilch),
		cbind(Zilch,A)
	),name="dV_dva2"),
	mxAlgebra(rbind(
		cbind(Zilch,Zilch),
		cbind(Zilch,Zilch)
	),name="dV_dmu"),
	mxAlgebra(rbind(
		cbind(Zilch,V11),
		cbind(V11,V11%x%(2*B))
	),name="dV_db"),
	mxAlgebra(rbind(Uno,Zip)+rbind(Zip,Uno%x%B),name="dyhat_dmu1"),
	mxAlgebra(rbind(Zip,Uno),name="dyhat_dmu2"),
	mxAlgebra(rbind(Zip,Zip),name="dyhat_dV"),
	
	mxData(observed=widedata,type="raw",sort=F),
	mxExpectationGREML(V="V",yvars=c("y1","y2"),REML=F,yhat="yhat"),
	mxFitFunctionGREML(
		dV=c(ve1="dV_dve1",ve2="dV_dve2",va1="dV_dva1",va2="dV_dva2",mu1="dV_dmu",mu2="dV_dmu",b="dV_db"),
		dyhat=c(ve1="dyhat_dV",ve2="dyhat_dV",va1="dyhat_dV",va2="dyhat_dV",mu1="dyhat_dmu1",mu2="dyhat_dmu2",b="M1")
	)
)
m3run <- mxRun(m3)
summary(m3run)

omxCheckCloseEnough(m1run$output$fit,m3run$output$fit,5e-5)
omxCheckCloseEnough(coef(m1run),coef(m3run)[1:4],0.01)
omxCheckCloseEnough(m1sum$GREMLfixeff$coeff,coef(m3run)[5:7],5e-4)
omxCheckWithinPercentError(m1run$output$standardErrors[1:4],m3run$output$standardErrors[1:4])
omxCheckCloseEnough(m1sum$GREMLfixeff$se,m3run$output$standardErrors[5:7],5e-5)
omxCheckCloseEnough(m1run$output$gradient,m2run$output$gradient[1:4],1e-5)




options(mxCondenseMatrixSlots=TRUE)
mxOption(reset=TRUE)
