require(OpenMx)

set.seed(1234)
dat <- cbind(rnorm(5),rep(1,5))
colnames(dat) <- c("y","x")
dat[5,1] <- NA

ge <- mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=FALSE)
gff <- mxFitFunctionGREML(dV=c(ve="I"))
plan <- mxComputeSequence(freeSet=c("Ve"),steps=list(
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
	mxMatrix(type="Unit", nrow=5,ncol=1,name="Uni"),
	mxAlgebra(vec2diag(Uni), name="I"),
	mxAlgebra( vec2diag(Uni%x%Ve), name="V"),
	ge,
	gff,
	plan
)

testrun <- mxRun(testmod)

summary(testrun)
omxCheckEquals(mxEval(nrow(I),testrun,T), 5)
omxCheckEquals(mxEval(ncol(I),testrun,T), 5)


testmod2 <- mxModel(
	"GREMLtest",
	mxData(observed = dat, type="raw", sort=FALSE),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
	mxMatrix(type="Unit", nrow=5,ncol=1,name="Uni"),
	mxAlgebra(vec2diag(Uni), name="I"),
	mxAlgebra( I%x%Ve, name="V"),
	ge,
	gff,
	plan
)

testrun2 <- mxRun(testmod2)

summary(testrun2)
omxCheckEquals(mxEval(nrow(I),testrun2,T), 5)
omxCheckEquals(mxEval(ncol(I),testrun2,T), 5)

