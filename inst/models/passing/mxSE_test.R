require(OpenMx)

set.seed(1611150)
x <- matrix(rnorm(1000,sd=2))
colnames(x) <- "x"

varmod <- mxModel(
	"mod",
	mxData(observed=x,type="raw"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0,labels="mu",name="Mu"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=4,labels="sigma2",name="Sigma2",lbound=0),
	mxExpectationNormal(covariance="Sigma2",means="Mu",dimnames=c("x")),
	mxAlgebra(sqrt(Sigma2),name="Sigma"),
	mxFitFunctionML()
)
varrun <- mxRun(varmod)

sdmod <- mxModel(
	"mod",
	mxData(observed=x,type="raw"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0,labels="mu",name="Mu"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=2,labels="sigma",name="Sigma",lbound=0),
	mxExpectationNormal(covariance="Sigma2",means="Mu",dimnames=c("x")),
	mxAlgebra(Sigma^2,name="Sigma2"),
	mxFitFunctionML()
)
sdrun <- mxRun(sdmod)

omxCheckCloseEnough(
	varrun$output$standardErrors[2],
	mxSE(x=Sigma2,model=sdrun),
	1e-6
)

omxCheckCloseEnough(
	sdrun$output$standardErrors[2],
	mxSE(x=Sigma,model=varrun),
	1e-6
)

