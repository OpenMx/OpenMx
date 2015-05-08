# Contributed by Charles Driver <driver at mpib-berlin.mpg.de>
#
# Based on:
#
# Manuel C. Voelkle, Johan H. L. Oud, Continuous time modelling with
# individually varying time intervals for oscillating and
# non-oscillating processes British Journal of Mathematical and
# Statistical Psychology Volume 66, Issue 1, pages 103-126.

library(OpenMx)
library(mvtnorm)

#mxOption(NULL, "Default optimizer", "NPSOL")

if (0) {
  set.seed(1)
  
  ############create data
  bivprocesscov<-matrix(c(
    0.63,0.24,0.45,0.24,0.33,0.23,0.2,0.2,0.14,0.17,
    0.24,0.46,0.24,0.42,0.23,0.38,0.2,0.32,0.17,0.26,
    0.45,0.24,0.65,0.24,0.46,0.24,0.26,0.22,0.17,0.19,
    0.24,0.42,0.24,0.52,0.25,0.47,0.23,0.39,0.2,0.33,
    0.33,0.23,0.46,0.25,0.65,0.25,0.34,0.24,0.21,0.21,
    0.23,0.38,0.24,0.47,0.25,0.57,0.26,0.47,0.24,0.39,
    0.2,0.2,0.26,0.23,0.34,0.26,0.67,0.28,0.36,0.26,
    0.2,0.32,0.22,0.39,0.24,0.47,0.28,0.64,0.29,0.53,
    0.14,0.17,0.17,0.2,0.21,0.24,0.36,0.29,0.68,0.31,
    0.17,0.26,0.19,0.33,0.21,0.39,0.26,0.53,0.31,0.69),nrow=10)
  
  bivprocessmeans<-c(2.5,2.84,2.57,2.84,2.6,2.84,2.65,2.84,2.66,2.84)
  
  subjects<-500
  bivprocess<-mvtnorm::rmvnorm(n=subjects,mean=bivprocessmeans,sigma=bivprocesscov)
  colnames(bivprocess)<-paste("V",1:10, sep="")
  
  intervals<-matrix(c(1,1,2,2),byrow=T,nrow=subjects,ncol=4)
  intervals<-matrix(rnorm(length(intervals),intervals,.1),ncol=4)
  
  colnames(intervals)<-paste("I",1:4, sep="")
  
  data <- round(cbind(bivprocess,intervals), 3)
  write.table(data, file="ct.csv", row.names=FALSE)
} else {
  data <- suppressWarnings(try(read.table("models/nightly/data/ContinuousTime.csv", header=TRUE), silent=TRUE))
  if (is(data, "try-error")) data <- read.table("data/ContinuousTime.csv", header=TRUE)
}


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
  EXPalgs[i]<-eval(substitute(mxAlgebra(theExpression,name = algName),
    list(theExpression = parse(text=fullAlgString)[[1]])))
}

### intercept (INT)
INTalgs<-list()
for(i in 1:4){
  algName<-paste("intd", i, sep="")
  defcall<-paste("data.I",i, sep="")
  fullAlgString<-paste("solve(DRIFT)%*%(omxExponential(DRIFT %x%",defcall,")-II)%*%t(CINT)", sep="")
  INTalgs[i]<-eval(substitute(mxAlgebra(theExpression,name = algName),
    list(theExpression = parse(text=fullAlgString)[[1]])))
}

### error covariance (Q)  
Qdalgs<-list()
for(i in 1:4){
  defcall<-paste("data.I",i, sep="")
  algName<-paste("Qd", i, sep="")
  fullAlgString<-paste("solve(DRIFTHATCH)%*%((omxExponential(DRIFTHATCH %x%",defcall,"))-(II%x%II))%*%rvectorize(Q)", sep="")
  Qdalgs[i]<-eval(substitute(mxAlgebra(theExpression,name = algName),list(theExpression = parse(text=fullAlgString)[[1]])))
}



#### OPENMX MODEL


larModel <-mxModel("lar",
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
#larModel <- mxOption(larModel, "Number of Threads", 1L)

larFit <- mxModel(larModel, mxComputeOnce('fitfunction', 'fit'))
got <- mxRun(larFit, silent=TRUE)
omxCheckCloseEnough(got$output$fit, 43550.90, .01)

larModel <- mxOption(larModel, "Feasibility tolerance", 2e-2)
testfit<-mxRun(larModel, intervals = TRUE)

omxCheckCloseEnough(testfit$output$fit, 8383.78, .1)
#cat(deparse(round(testfit$output$estimate,3)))
est <- c(0.675, 0.266, 0.466, 2.496, 2.858, -0.464, 0.056,  0.246, -0.131, 0.526, -0.014, 0.164, 0.524, 0.224)
if (0) {
  testfit$output$estimate - est
  max(abs(testfit$output$estimate - est))
}
omxCheckCloseEnough(testfit$output$estimate, est, .01)

ci <- testfit$output$confidenceIntervals
#cat(deparse(round(ci[,'ubound'], 4)))
omxCheckCloseEnough(ci[,'lbound'], c(-0.518, 0.1997, 0.0282, -0.1568), .03)
omxCheckCloseEnough(ci[,'ubound'], c(-0.4135, 0.2948, 0.0837, -0.1062), .03)

if (0) {
  omxCheckCloseEnough(testfit$output$iterations, 9, 1)
  omxCheckCloseEnough(testfit$output$status$code, -1)
}
