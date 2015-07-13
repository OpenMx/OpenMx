require(OpenMx)

omxCheckError(mxExpectationGREML(V=1),
              "argument 'V' is not of type 'character' (the name of the expected covariance matrix)")

set.seed(1234)
dat <- cbind(rnorm(100),rep(1,100))
colnames(dat) <- c("y","x")
testmod <- mxModel(
  "GREMLtest",
  #mxData(observed=dat, type="raw", sort=F),
  mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
  mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
  mxAlgebra(I %x% Ve,name="V"),
  mxExpectationGREML(V="V",Xvars=list("x"),yvars="y",addOnes=F),
  mxFitFunctionGREML()
)
omxCheckError(mxRun(testmod),
              "the GREML expectation function does not have a dataset associated with it in model 'GREMLtest'")

testmod <- mxModel(
  testmod,
  mxData(observed=dat, type="raw", sort=F),
  mxMatrix("Full",1,1,F,labels="data.y",name="Z")
)
omxCheckError(mxRun(testmod),
              "definition variables are incompatible (and unnecessary) with GREML expectation")

testmod <- mxModel(
  "GREMLtest",
  mxData(observed=dat, type="raw", sort=F),
  mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
  mxMatrix("Full",nrow=100,ncol=99,free=F,values=diag(100)[,-100],name="V",condenseSlots=T),
  mxExpectationGREML(V="V",dataset.is.yX=T),
  mxFitFunctionGREML()
)
omxCheckError(mxRun(testmod),
              "'V' matrix is not square")

testmod$V <- mxMatrix("Full",nrow=99,ncol=99,free=F,values=diag(100)[-100,-100],name="V",condenseSlots=T)
omxCheckError(mxRun(testmod),
              "y and V matrices do not have equal numbers of rows")

testmod <- mxModel(
  "GREMLtest",
  mxData(observed=dat, type="raw", sort=F),
  mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 1, labels = "ve", lbound = 0.0001, name = "Ve"),
  mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
  mxAlgebra(I %x% Ve,name="V"),
  mxExpectationGREML(V="V",dataset.is.yX=TRUE,casesToDropFromV=101L),
  mxFitFunctionGREML()
)
omxCheckWarning(
  mxRun(testmod, suppressWarnings=TRUE),
  "casesToDrop vector in GREML expectation contains indices greater than the number of datapoints")

testmod <- mxModel(
  "GREMLtest",
  mxData(observed=dat, type="raw", sort=F),
  mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
  mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
  mxMatrix("Iden",nrow=99,name="J",condenseSlots=T),
  mxAlgebra(I %x% Ve,name="V"),
  mxExpectationGREML(V="V", dataset.is.yX = T),
  mxFitFunctionGREML(dV = c(ve="J"))
)
omxCheckError(mxRun(testmod),
              "all derivatives of V must have the same dimensions as V")

testmod <- mxModel(
  "GREMLtest",
  mxData(observed=dat, type="raw", sort=F),
  mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
  mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
  mxAlgebra(I %x% Ve,name="V"),
  mxExpectationGREML(V="V", dataset.is.yX=T),
  mxFitFunctionGREML(dV = c(ve="J"))
)
omxCheckError(mxRun(testmod),
              "The reference 'J' does not exist.  It is used by named reference 'GREMLtest.fitfunction' .")


testmod <- mxModel(
  "GREMLtest",
  mxData(observed = matrix(dat[,1],1,100,dimnames=list(NULL,paste("y",1:100,sep=""))), type="raw"),
  mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
  mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
  mxMatrix("Zero",1,100,name="Zm"),
  mxAlgebra(I %x% Ve,name="V"),
  mxExpectationNormal(covariance="V",means="Zm", dimnames=paste("y",1:100,sep="")),
  mxFitFunctionGREML()
)
omxCheckError(mxRun(testmod),
              "GREML fitfunction is currently only compatible with GREML expectation")


testmod <- mxModel(
	"GREMLtest",
	mxData(observed=dat, type="raw", sort=F),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
	mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
	mxAlgebra(I %x% Ve,name="V"),
	mxExpectationGREML(V="V",Xvars=list("x"),yvars="y",addOnes=F),
	mxFitFunctionGREML()
)
omxCheckError(mxRefModels(testmod),
							"reference models for GREML expectation not implemented")
testmod <- mxModel(
	"GREMLtest",
	mxData(observed=dat, type="raw", sort=F),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
	mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
	mxAlgebra(I %x% Ve,name="V"),
	mxExpectationGREML(V="V",Xvars=list("x"),yvars="y",addOnes=F),
	mxFitFunctionML()
)
omxCheckError(mxRefModels(testmod),
							"reference models for GREML expectation not implemented")


testmod <- mxModel(
	"GREMLtest",
	mxData(observed=dat, type="raw", sort=F),
	mxMatrix(type = "Unit", nrow = 100, ncol=100, name = "V", condenseSlots = T),
	mxExpectationGREML(V="V",Xvars=list("x"),yvars="y",addOnes=F),
	mxFitFunctionGREML()
)
omxCheckError(mxRun(testmod),
							"Expected covariance matrix is non-positive-definite at initial values")


z <- matrix(-1,100,2)
colnames(z) <- c("z1","z2")
dat2 <- cbind(dat,z)
testmod <- mxModel(
	"GREMLtest",
	mxData(observed=dat2, type="raw", sort=F),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
	mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
	mxAlgebra(I %x% Ve,name="V"),
	mxExpectationGREML(V="V",Xvars=list(c("x","z1","z2")),yvars="y",addOnes=F),
	mxFitFunctionGREML()
)
omxCheckError(mxRun(testmod),
	"Cholesky factorization failed at initial values; possibly, the matrix of covariates is rank-deficient")

