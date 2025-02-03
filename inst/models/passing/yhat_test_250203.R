require(OpenMx)

set.seed(1234)
dat <- cbind(rnorm(100),rep(1,100))
colnames(dat) <- c("y","x")

ge <- mxExpectationGREML(V="V",yvars="y",REML=FALSE,yhat="foo")
gff <- mxFitFunctionGREML(dV=c(ve="I"))
plan <- mxComputeSequence(steps=list(
	mxComputeNewtonRaphson(fitfunction="fitfunction"),
	mxComputeOnce('fitfunction', c('fit','gradient','hessian','ihessian')),
	mxComputeStandardError(),
	mxComputeReportDeriv(),
	mxComputeReportExpectation()
))

testmod <- mxModel(
	"GREMLtest",
	mxData(observed = dat, type="raw", sort=FALSE),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
	mxMatrix(type="Full",nrow=100,ncol=1,name="foo",condenseSlots=T,free=F,values=0.12345),
	mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
	mxAlgebra(I %x% Ve,name="V"),
	ge,
	gff,
	plan
)

omxCheckError(
	mxRun(testmod),
	"the first element of the yhat vector is 0.123450"
)
