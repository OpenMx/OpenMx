require(OpenMx)

if (mxOption(NULL,"Default optimizer")!='NPSOL') stop("SKIP")

set.seed(1611150)
x <- matrix(rnorm(1000))
colnames(x) <- "x"

m1 <- mxModel(
	"m1",
	mxData(observed=x,type="raw"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0,labels="mu",name="Mu"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=1,labels="sigma2",name="Sigma",lbound=0),
	mxExpectationNormal(covariance="Sigma",means="Mu",dimnames=c("x")),
	mxFitFunctionML()
)
run1 <- mxRun(m1)

hess.ini <- matrix(-2*c(-1000/1,0,0,-1000/2),nrow=2,ncol=2)
ws <- chol(hess.ini)

plan <- mxComputeSequence(steps=list(
	mxComputeGradientDescent(engine="NPSOL",warmStart=ws),
	mxComputeNumericDeriv(),
	mxComputeStandardError(),
	mxComputeHessianQuality(),
	mxComputeReportDeriv(),
	mxComputeReportExpectation()
))

m2 <- mxModel(
	"m2",
	plan,
	mxData(observed=x,type="raw"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0,labels="mu",name="Mu"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=1,labels="sigma2",name="Sigma",lbound=0),
	mxExpectationNormal(covariance="Sigma",means="Mu",dimnames=c("x")),
	mxFitFunctionML()
)
run2 <- mxRun(m2)

if(run2$output$evaluations > run1$output$evaluations){stop("WarmStartTest.R failed")}

plan2 <- mxComputeSequence(steps=list(
	mxComputeGradientDescent(engine="NPSOL",warmStart=matrix(1,1,1)),
	mxComputeNumericDeriv(),
	mxComputeStandardError(),
	mxComputeHessianQuality(),
	mxComputeReportDeriv(),
	mxComputeReportExpectation()
))

m3 <- mxModel(
	"m3",
	plan2,
	mxData(observed=x,type="raw"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=-0.3,labels="mu",name="Mu"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0.97,labels="sigma2",name="Sigma",lbound=0),
	mxExpectationNormal(covariance="Sigma",means="Mu",dimnames=c("x")),
	mxFitFunctionML()
)
omxCheckWarning(mxRun(m3),"warmStart size 1 does not match number of free parameters 2 (ignored)")

