library(OpenMx)
if(mxOption(NULL,"Default optimizer")!="SLSQP"){stop("SKIP")}

mxOption(NULL,"Nudge zero starts","No")
mxOption(NULL,"Calculate Hessian","No")
mxOption(NULL,"Standard Errors","No")

mxd <- data.frame(
	y=mxFactor(as.numeric(runif(100)<.5), levels=c(0,1)),
	k5=runif(100),
	k618=runif(100),
	age=round(rnorm(100,30,5)),
	wc=rbinom(100,1,.5),
	hc=rbinom(100,1,.5),
	lwg=rnorm(100))

mxd2 <- mxd
mxd2$y <- as.numeric(mxd2$y)-1

mod2b <- mxModel(
	"constrLPM",
	mxData(observed=mxd2,type="raw"),
	mxMatrix(
		type="Full",nrow=1,ncol=7,free=F,labels=c("data.k5","data.k618","data.age","data.wc","data.hc","data.lwg",NA),name="X"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,
					 values=0,
					 labels="a",name="A"),
	mxMatrix(type="Full",nrow=7,ncol=1,free=T,
					 values=0,
					 labels=paste("b",1:7,sep=""),name="B"),
	mxAlgebra(A + X%*%B, name="Yhat"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=1,labels="sigma",lbound=0.0001,name="Sigma"),
	mxExpectationNormal(covariance="Sigma",means="Yhat",dimnames=c("y")),
	mxFitFunctionML(),
	mxConstraint(-1%x%Yhat < 0, name="c1"),
	mxConstraint(Yhat - 1 < 0, name="c2")
)
mod2b <- omxCheckWarning(
	mxRun(mod2b),
  "Constraint 'constrLPM.c2' depends on definition variables; This may not do what you expect. See ?mxConstraint")
