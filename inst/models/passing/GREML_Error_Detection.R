#
#   Copyright 2007-2018 The OpenMx Project
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


require(OpenMx)
options(mxCondenseMatrixSlots=TRUE)
mxOption(NULL,"Analytic Gradients","Yes")
require(mvtnorm)

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
							"Reference models for GREML expectation are not implemented")



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


set.seed(476)
A1 <- matrix(0,100,100)  
A1[lower.tri(A1)] <- runif(4950, -0.025, 0.025)
A1 <- A1 + t(A1)
diag(A1) <- runif(100,0.95,1.05)
A2 <- matrix(0,100,100)  
A2[lower.tri(A2)] <- runif(4950, -0.025, 0.025)
A2 <- A2 + t(A2)
diag(A2) <- runif(100,0.95,1.05)
y <- t(rmvnorm(1,sigma=A1*0.25)+rmvnorm(1,sigma=A2*0.25))  
y <- y + rnorm(100,sd=sqrt(0.5))
x <- rnorm(100) 
dat3 <- cbind(y,x)
rm(x,y)
colnames(dat3) <- c("y","x")
testmod <- mxModel(
	"GREMLtest",
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values =0.5, labels = "ve", lbound = 0.0001, 
					 name = "Ve"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va1", name = "Va1"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va2", name = "Va2"),
	mxData(observed = dat3, type="raw", sort=FALSE),
	mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T),
	mxMatrix("Iden",nrow=100,name="I"),
	mxMatrix("Symm",nrow=100,free=F,values=A1,name="A1"),
	mxMatrix("Symm",nrow=100,free=F,values=A2,name="A2"),
	mxAlgebra((A1%x%Va1) + (A2%x%Va2) + (I%x%Ve), name="V"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=F,values=0.64,name="aug"),
	mxMatrix(type="Zero",nrow=1,ncol=1,name="Zilch"),
	mxFitFunctionGREML(dV=c(ve="I",va1="A1",va2="A2"),aug="aug",augHess="Zilch")
)
omxCheckError(
	mxRun(testmod),
	"if argument 'augHess' has nonzero length, then argument 'augGrad' must as well")


testmod$fitfunction <- mxFitFunctionGREML(dV=c(ve="I",va1="A1",va2="A2"),aug="aug")
omxCheckError(
	mxRun(testmod),
	"if arguments 'dV' and 'aug' have nonzero length, then 'augGrad' must as well")


testmod$fitfunction <- mxFitFunctionGREML(dV=c(ve="I",va1="A1",va2="A2",va3="I"))
omxCheckError(
	mxRun(testmod),
	"Problem in dVnames mapping")


testmod <- mxModel(
	"GREMLtest",
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values =0.5, labels = "ve", lbound = 0.0001, 
					 name = "Ve"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va1", name = "Va1"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va2", name = "Va2"),
	mxData(observed = dat3, type="raw", sort=FALSE),
	mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T),
	mxMatrix("Iden",nrow=100,name="I"),
	mxMatrix("Symm",nrow=100,free=F,values=A1,name="A1"),
	mxMatrix("Symm",nrow=100,free=F,values=A2,name="A2"),
	mxAlgebra((A1%x%Va1) + (A2%x%Va2) + (I%x%Ve), name="V"),
	mxComputeSequence(steps=list(
		mxComputeNewtonRaphson(fitfunction="fitfunction"),
		mxComputeOnce('fitfunction', c('fit','gradient','hessian','ihessian')),
		mxComputeStandardError(),
		mxComputeReportDeriv(),
		mxComputeReportExpectation()
	)),
	mxFitFunctionGREML(dV=c(ve="I",va1="A1"))
)
omxCheckError(
	mxRun(testmod),
	"At least one free parameter has no corresponding element in 'dV'")


testmod$compute <- mxComputeDefault()
omxCheckError(
	mxRun(testmod),
	"At least one free parameter has no corresponding element in 'dV'")



testmod <- mxModel(
	"GREMLtest",
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values =0.5, labels = "ve", lbound = 0.0001, 
					 name = "Ve"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va1", name = "Va1"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va2", name = "Va2"),
	mxData(observed = dat3, type="raw", sort=FALSE),
	mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T),
	mxMatrix("Iden",nrow=100,name="I"),
	mxMatrix("Symm",nrow=100,free=F,values=A1,name="A1"),
	mxMatrix("Symm",nrow=100,free=F,values=A2,name="A2"),
	mxAlgebra((A1%x%Va1) + (A2%x%Va2) + (I%x%Ve), name="V"),
	mxComputeSequence(steps=list(
		mxComputeNewtonRaphson(fitfunction="fitfunction"),
		mxComputeOnce('fitfunction', c('fit','gradient','hessian','ihessian')),
		mxComputeStandardError(),
		mxComputeReportDeriv(),
		mxComputeReportExpectation()
	)),
	mxFitFunctionGREML(dV=c(ve="I",va1="A1",va2="A2",va3="I"))
)
omxCheckError(
	mxRun(testmod),
	"Problem in dVnames mapping")



testmod <- mxModel(
	"GREMLtest",
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values =0.5, labels = "ve", lbound = 0.0001, 
					 name = "Ve"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va1", name = "Va1"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va2", name = "Va2"),
	mxData(observed = dat3, type="raw", sort=FALSE),
	mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T),
	mxMatrix("Iden",nrow=100,name="I"),
	mxMatrix("Symm",nrow=100,free=F,values=A1,name="A1"),
	mxMatrix("Symm",nrow=100,free=F,values=A2,name="A2"),
	mxAlgebra((A1%x%Va1) + (A2%x%Va2) + (I%x%Ve), name="V"),
	mxComputeSequence(steps=list(
		mxComputeNewtonRaphson(fitfunction="fitfunction"),
		mxComputeOnce('fitfunction', c('fit','gradient','hessian','ihessian')),
		mxComputeStandardError(),
		mxComputeReportDeriv(),
		mxComputeReportExpectation()
	)),
	mxMatrix(type="Full",nrow=1,ncol=1,free=F,values=2.1,name="aug"),
	mxMatrix(type="Zero",nrow=2,ncol=1,free=F,name="ag"),
	mxMatrix(type="Zero",nrow=3,ncol=3,free=F,name="ah"),
	mxFitFunctionGREML(dV=c(ve="I",va1="A1",va2="A2"),aug="aug",augGrad="ag",augHess="ah")
)
omxCheckError(
	mxRun(testmod),
	"matrix referenced by 'augGrad' must have same number of elements as argument 'dV'")



testmod$ag <- mxMatrix(type="Zero",nrow=3,ncol=1,free=F,name="ag")
testmod$ah <- mxMatrix(type="Zero",nrow=2,ncol=3,free=F,name="ah")
omxCheckError(
	mxRun(testmod),
	"matrix referenced by 'augHess' must be square (instead of 2x3)")



testmod$ah <- mxMatrix(type="Zero",nrow=2,ncol=2,free=F,name="ah")
omxCheckError(
	mxRun(testmod),
	"Augmentation derivatives non-conformable (gradient is size 3 and Hessian is 2x2)")

