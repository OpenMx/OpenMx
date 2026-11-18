#
#   Copyright 2007-2026 by the individuals mentioned in the source code history
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

#Based on models/nightly/ContinuousTime.R by Charles Driver

library(OpenMx)
#Use only CSOLNP for this test, because NPSOL is incompatible with parallelized CIs, 
#and SLSQP can parallelize its numeric gradient elements if it is prevented from parallelizing CIs or rowwise fit evaluations:
if(mxOption(NULL,"Default optimizer")!="CSOLNP"){stop("SKIP")}
#mxOption(key='Number of Threads', value=4)
#mxOption(NULL, "Default optimizer", "CSOLNP")
ompThreads = suppressWarnings(as.integer(Sys.getenv("OMP_NUM_THREADS", unset = "4")))
if (is.na(ompThreads) || ompThreads < 1) ompThreads = 4L
nthreads = min(ompThreads, 4L)
if (nthreads < 2) stop("SKIP")
library(mvtnorm)


data <- suppressWarnings(try(read.table("models/nightly/data/ContinuousTime.csv", header=TRUE), silent=TRUE))
if (is(data, "try-error")) data <- read.table("data/ContinuousTime.csv", header=TRUE)


#########create variables and algebras for openmx

### VARIABLES
manifestV1<-paste("V",seq(1,2*5,2), sep="")
manifestV2<-paste("V",seq(2,2*5,2), sep="")

latentV1<-paste("L1_",seq(1,5,1), sep="")
latentV1r<-paste(latentV1[-1],"r", sep="")

latentV2<-paste("L2_",seq(1,5,1), sep="")
latentV2r<-paste(latentV2[-1],"r", sep="")


### ALGEBRAS

###DRIFT
EXPalgs<-list()
for(i in 1:4){
	defcall<-paste("data.I",i, sep="")
	algName<-paste("EXPd", i, sep="")
	fullAlgString<-paste("omxExponential(DRIFT %x%",defcall,")", sep="")
	EXPalgs[i]<- list(eval(substitute(mxAlgebra(theExpression,name = algName),
																		list(theExpression = parse(text=fullAlgString)[[1]]))))
}

### intercept (INT)
INTalgs<-list()
for(i in 1:4){
	algName<-paste("intd", i, sep="")
	defcall<-paste("data.I",i, sep="")
	fullAlgString<-paste("solve(DRIFT)%*%(omxExponential(DRIFT %x%",defcall,")-II)%*%t(CINT)", sep="")
	INTalgs[i]<-list(eval(substitute(mxAlgebra(theExpression,name = algName),
																	 list(theExpression = parse(text=fullAlgString)[[1]]))))
}

### error covariance (Q)  
Qdalgs<-list()
for(i in 1:4){
	defcall<-paste("data.I",i, sep="")
	algName<-paste("Qd", i, sep="")
	fullAlgString<-paste("solve(DRIFTHATCH)%*%((omxExponential(DRIFTHATCH %x%",defcall,"))-(II%x%II))%*%rvectorize(Q)", sep="")
	Qdalgs[i]<-list(eval(substitute(mxAlgebra(theExpression,name = algName),list(theExpression = parse(text=fullAlgString)[[1]]))))
}



#### OPENMX MODEL


larModel <-mxModel(
	"lar",
	mxData(data,type="raw"),  ####################### Modify this line to change data file
	manifestVars=c(manifestV1,manifestV2),
	latentVars=c(latentV1,"i1","L1t",latentV1r,
							 latentV2,"i2","L2t",latentV2r),
	
	EXPalgs,
	Qdalgs,
	INTalgs,
	
	mxMatrix(type="Iden", nrow=2, ncol=2, free=FALSE, name="II"),
	mxMatrix("Full", values=(matrix(1,2,2)-diag(2)), name = "tempa"),
	
	mxMatrix(type="Full", labels=c("F11","F12","F21","F22"), values=c(-.1,.01,.01,-.2), free=T,nrow=2,ncol=2, name="DRIFT"),
	mxCI("DRIFT"),
	
	mxMatrix(type="Full", labels=c("Q11","Q21","Q21","Q22"), values=c(.5,.1,.1,.5), free=T, name="Q",nrow=2,ncol=2),    
	mxAlgebra(DRIFT%x%II + II%x%DRIFT, name = "DRIFTHATCH"),
	mxMatrix("Full", values=(matrix(1,4,4)-diag(4)), name = "tempb"),   
	
	mxMatrix(type="Full", labels=c("cint1","cint2"), values=c(2,2), free=T, name="CINT",nrow=1,ncol=2),
	
	#latent to manifest loadings
	mxPath(from=c(latentV1,latentV2),
				 to=c(manifestV1,manifestV2),
				 connect="single",arrows=1,free=F,values=1),
	
	#manifest means
	mxPath(from="one",   #set manifest means to 0
				 to=c(manifestV1,manifestV2),
				 arrows=1,free=FALSE,values=0),
	
	#   add CINT   
	mxPath(to=c(latentV1[-1]),
				 from=c("L1t"),
				 values=NA,free=F,arrows=1,connect="single",labels=paste("intd",1:4,"[",1,",1]", sep="")),
	mxPath(to=c(latentV2[-1]),
				 from=c("L2t"),
				 values=NA,free=F,arrows=1,connect="single",labels=paste("intd",1:4,"[",2,",1]", sep="")),
	
	# continuous intercept mean
	mxPath(from="one",to=c("L1t","L2t"),arrows=1,free=F,values=1),    
	
	#initial intercept loadings
	mxPath(from="i1",
				 to=latentV1[1],
				 arrows=1,values=1,free=F),
	
	mxPath(from="i2",
				 to=latentV2[1],
				 arrows=1,values=1,free=F),
	
	#initial intercept means
	mxPath(from="one",
				 to=c("i1","i2"),
				 arrows=1,free=c(T,T),values=c(.5,.5),
				 labels=c("m1","m2")),
	
	#intercept variance and cov (phi)
	mxPath(from=c("i1","i2"),
				 connect="unique.pairs",arrows=2,free=T,values=c(.5,.1,.5),
				 labels=c("phi11","phi21","phi22")),
	
	#AR
	mxPath(from=latentV1[-length(latentV1)],
				 to=latentV1[-1],
				 arrows=1,free=F,values=NA,connect="single",labels=paste("EXPd",1:4,"[1,1]", sep="")),
	mxPath(from=latentV2[-length(latentV2)],
				 to=latentV2[-1],
				 arrows=1,free=F,values=NA,connect="single",labels=paste("EXPd",1:4,"[2,2]", sep="")),
	
	# # CL 
	mxPath(from=latentV1[-length(latentV2)],to=latentV2[-1],
				 arrows=1,free=F,values=NA, labels=paste("EXPd",1:4,"[2,1]", sep="")),
	
	mxPath(from=latentV2[-length(latentV2)],to=latentV1[-1],
				 arrows=1,free=F,values=NA,labels=paste("EXPd",1:4,"[1,2]", sep="")),
	
	#add latent disturbances
	mxPath(to=latentV1[-1],from=latentV1r,
				 connect="single",free=F,arrows=1,values=1),
	mxPath(to=latentV2[-1],from=latentV2r,
				 connect="single",arrows=1,free=F,values=1),
	
	#latent disturbance variance
	mxPath(from=latentV1r,free=F,values=NA,arrows=2,labels=paste("Qd",1:4,"[1,1]", sep="")),
	mxPath(from=latentV2r,arrows=2,free=F,values=NA,labels=paste("Qd",1:4,"[4,1]", sep="")),
	
	#disturbance cross correlations between same time point
	mxPath(from=c(latentV1r,latentV2r),to=c(latentV2r,latentV1r),
				 arrows=2,free=F,values=NA,connect="single",labels=paste("Qd",1:4,"[2,1]", sep="")),
	type="RAM") 

larModel <- mxOption(larModel, "Calculate Hessian", "No")
larModel <- mxOption(larModel, "Standard Errors", "No")
larModel$fitfunction$rowwiseParallel <- FALSE

larFit <- mxModel(larModel, mxComputeOnce('fitfunction', 'fit'))
got <- mxRun(larFit, silent=TRUE)
omxCheckCloseEnough(got$output$fit, 43550.90, .01)

larModel <- mxOption(larModel, "Feasibility tolerance", 2e-2)
mxOption(NULL,"Parallelize confidence intervals","Yes")
testfit1 <- mxRun(larModel, intervals = TRUE)
mxOption(NULL,"Parallelize confidence intervals","No")
testfit2 <- mxRun(larModel, intervals = TRUE)

omxCheckCloseEnough(testfit1$output$confidenceIntervals,testfit2$output$confidenceIntervals,1e-7)
( summ1 <- summary(testfit1) )
( summ2 <- summary(testfit2) )
omxCheckTrue(summ1$wallTime < summ2$wallTime)




mxOption(reset=TRUE)
