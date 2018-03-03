library(OpenMx)
if(mxOption(NULL,"Default optimizer")!="SLSQP"){stop("SKIP")}

mxOption(NULL,"Nudge zero starts","No")
library(car)

head(Mroz)

t <- ifelse(Mroz[,1]=="yes",1,0)

mxd <- data.frame(
	y=mxFactor(t,levels=c(0,1)),
	k5=as.numeric(Mroz$k5),
	k618=as.numeric(Mroz$k618),
	age=as.numeric(Mroz$age),
	wc=as.numeric(Mroz$wc)-1,
	hc=as.numeric(Mroz$hc)-1,
	lwg=Mroz$lwg,inc=Mroz$inc)

mxd2 <- mxd
mxd2$y <- as.numeric(mxd2$y)-1

mod2b <- mxModel(
	"constrLPM",
	mxData(observed=mxd2,type="raw"),
	mxMatrix(
		type="Full",nrow=1,ncol=7,free=F,labels=c("data.k5","data.k618","data.age","data.wc","data.hc","data.lwg","data.inc"),name="X"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,
					 values=mean(t),
					 labels="a",name="A"),
	mxMatrix(type="Full",nrow=7,ncol=1,free=T,
					 values=0,
					 labels=paste("b",1:7,sep=""),name="B"),
	mxAlgebra(A + X%*%B, name="Yhat"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=var(t),labels="sigma",lbound=0.0001,name="Sigma"),
	mxExpectationNormal(covariance="Sigma",means="Yhat",dimnames=c("y")),
	mxFitFunctionML(),
	mxConstraint(-1%x%Yhat < 0, name="c1"),
	mxConstraint(Yhat - 1 < 0, name="c2")
)
mod2b <- omxCheckWarning(mxRun(mod2b),
                         "Constraint 'constrLPM.c2' depends on definition variables; This may not do what you expect. See ?mxConstraint")
