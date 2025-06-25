#
#   Copyright 2007-2020 by the individuals mentioned in the source code history
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


library(OpenMx)
library(testthat)
context("GREML Error Detection")

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


omxCheckError(mxModel(
	"GREMLtest",
	mxData(observed=dat, type="raw", sort=F),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
	mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
	mxAlgebra(I %x% Ve,name="V"),
	mxExpectationGREML(V="V",Xvars=list("x"),yvars="y",addOnes=F),
	mxFitFunctionGREML(autoDerivType="muneric")
),"'muneric' should be one of 'semiAnalyt' and 'numeric'")


testmod <- mxModel(
	"GREMLtest",
	mxData(observed=dat, type="raw", sort=F),
	mxMatrix(type = "Unit", nrow = 100, ncol=100, name = "V", condenseSlots = T),
	mxExpectationGREML(V="V",Xvars=list("x"),yvars="y",addOnes=F),
	mxFitFunctionGREML()
)
omxCheckError(mxRun(testmod),
							"Expected covariance matrix is non-positive-definite at initial values")



testmod <- mxModel(
	"GREMLtest",
	mxData(observed=dat, type="raw", sort=F),
	mxMatrix(type="Diag", nrow=100, ncol=100, free=T, values=0.1, labels="v", lbound=0.0001, name="V"),
	mxExpectationGREML(V="V",Xvars=list("x"),yvars="y",addOnes=F),
	mxFitFunctionML()
)
omxCheckWarning(
	mxRun(testmod),
	"use of an ML fitfunction with a GREML expectation is deprecated; instead, try using a GREML fitfunction, with argument `REML=FALSE` provided to `mxExpectationGREML()`")



testmod <- mxModel(
	"GREMLtest",
	mxData(observed=dat, type="raw", sort=F),
	mxMatrix(type = "Unit", nrow = 100, ncol=100, name = "V", condenseSlots = T),
	mxExpectationGREML(V="V",Xvars=list("x"),yvars="y",addOnes=F),
	mxFitFunctionML()
)
omxCheckError(mxRun(testmod),
							"Expected covariance matrix is non-positive-definite at initial values")



testmod <- mxModel(
	"GREMLtest",
	mxData(observed = dat, type="raw", sort=FALSE),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
	mxMatrix(type="Full",nrow=100,ncol=1,name="foo",free=T,values=0.12345,labels="bar"),
	mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
	mxAlgebra(I %x% Ve,name="V"),
	mxExpectationGREML(V="V",yvars="y",REML=FALSE,yhat="foo"),
	mxFitFunctionGREML(autoDerivType="semiAnalyt")
)


testmod <- mxModel(
	"GREMLtest",
	mxData(observed = dat, type="raw", sort=FALSE),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
	mxMatrix(type="Full",nrow=100,ncol=1,name="foo",free=T,values=0.12345,labels="bar"),
	mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
	mxAlgebra(I %x% Ve,name="V"),
	mxExpectationGREML(V="V",yvars="y",REML=FALSE,yhat="foo"),
	mxFitFunctionGREML(dV=c(ve="I"),autoDerivType="numeric")
)


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


testmod <- mxModel(
	"GREMLtest",
	mxData(observed=dat2, type="raw", sort=F),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
	mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
	mxAlgebra(I %x% Ve,name="V"),
	mxExpectationGREML(V="V",Xvars=list(c("x","z1","z2")),yvars="y",addOnes=F),
	mxFitFunctionML()
)
omxCheckError(mxRun(testmod),
							"Cholesky factorization failed at initial values; possibly, the matrix of covariates is rank-deficient")


testmod <- mxModel(
	"GREMLtest",
	mxData(observed=dat, type="raw", sort=F),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
	mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
	mxAlgebra(I %x% Ve,name="V"),
	mxExpectationGREML(V="V",yvars="y",addOnes=F),
	mxFitFunctionGREML()
)
omxCheckWarning(
	mxRun(testmod),
	"argument 'addOnes' is FALSE, but no covariates are named in argument 'Xvars'; the 'X' matrix will be constructed for intercept(s)-only"
)

testmod <- mxModel(
	"GREMLtest",
	mxData(observed=dat, type="raw", sort=F),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
	mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
	mxAlgebra(I %x% Ve,name="V"),
	mxExpectationGREML(V="V",yvars="y",addOnes=F),
	mxFitFunctionGREML(dV=c(ve="I",ve="I"))
)
omxCheckError(
	mxRun(testmod),
	"duplicated element names in argument 'dV'"
)

testmod <- mxModel(
	"GREMLtest",
	mxData(observed=dat, type="raw", sort=F),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
	mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
	mxAlgebra(I %x% Ve,name="V"),
	mxExpectationGREML(V="V",yvars="y",Xvars=list("x"),addOnes=F),
	mxFitFunctionGREML(dV=c(ve="I"))
)
testmod@fitfunction@.parallelDerivScheme <- 4L
omxCheckWarning(
	mxRun(testmod),
	"`.parallelDerivScheme` not equal to 1, 2, or 3"
)


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
	"length of argument 'dV' is greater than the number of explicit free parameters")


testmod <- mxModel(
	"GREMLtest",
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values =0.5, labels = "ve", lbound = 0.0001, 
					 name = "Ve"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va1", name = "Va1"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va2", name = "Va2"),
	mxData(observed = dat3, type="raw", sort=FALSE),
	mxExpectationGREML(V="V",yvars="y",REML=F,yhat="foo"),
	mxMatrix("Iden",nrow=100,name="I"),
	mxMatrix("Symm",nrow=100,free=F,values=A1,name="A1"),
	mxMatrix("Symm",nrow=100,free=F,values=A2,name="A2"),
	mxAlgebra((A1%x%Va1) + (A2%x%Va2) + (I%x%Ve), name="V"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=F,values=0.64,name="aug"),
	mxMatrix(type="Zero",nrow=1,ncol=1,name="Zilch"),
	mxMatrix(type="Full",nrow=100,ncol=1,name="foo",free=T,values=0.12345,labels="bar"),
	mxMatrix(type="Unit",nrow=100,ncol=1,name="Uno"),
	mxFitFunctionGREML(dyhat=c(bar="Uno"),aug="aug",augHess="Zilch")
)
omxCheckError(
	mxRun(testmod),
	"if argument 'augHess' has nonzero length, then argument 'augGrad' must as well")

testmod$fitfunction <- mxFitFunctionGREML(dyhat=c(bar="Uno"),aug="aug")
omxCheckError(
	mxRun(testmod),
	"if arguments 'dyhat' and 'aug' have nonzero length, then 'augGrad' must as well"
)

testmod$fitfunction <- mxFitFunctionGREML(aug="aug")
omxCheckError(
	mxRun(testmod),
	"if using semi-analytic derivatives and 'aug' has nonzero length, then 'augGrad' must as well"
)


testmod <- mxModel(
	"GREMLtest",
	mxMatrix(
		type = "Full", nrow = 1, ncol=1, free=T, values =0.5, labels = "ve", lbound = 0.0001, 
		name = "Ve"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va1", name = "Va1"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va2", name = "Va2"),
	mxData(observed = dat3, type="raw", sort=FALSE),
	mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T),
	mxMatrix("Iden",nrow=100,name="I"),
	mxMatrix("Symm",nrow=100,free=F,values=A1,name="A1"),
	mxMatrix("Symm",nrow=100,free=F,values=A2,name="A2"),
	mxAlgebra((A1%x%Va1) + (A2%x%Va2) + (I%x%Ve), name="V"),
	mxFitFunctionGREML(dV=c(va1="A1",va2="A2"),autoDerivType="numeric")
)


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
	"length of argument 'dV' is greater than the number of explicit free parameters")



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
		mxComputeOnce('fitfunction', c('fit','gradient','hessian')),
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
	"matrix referenced by 'augGrad' must have as many elements as there are explicit free parameters")



testmod$ag <- mxMatrix(type="Zero",nrow=3,ncol=1,free=F,name="ag")
testmod$ah <- mxMatrix(type="Zero",nrow=2,ncol=3,free=F,name="ah")
omxCheckError(
	mxRun(testmod),
	"matrix referenced by 'augHess' must be square (instead of 2x3)")



testmod$ah <- mxMatrix(type="Zero",nrow=2,ncol=2,free=F,name="ah")
omxCheckError(
	mxRun(testmod),
	"Augmentation derivatives non-conformable (gradient is size 3 and Hessian is 2x2)")

####
# Test error messages relating to `mxExpectationGREML()`'s arguments, especially `REML`, especially `REML=FALSE`: ####
####

omxCheckError(
	mxExpectationGREML(V="V",yvars=10L),
	"argument 'yvars' is not of type 'character' (the data column names of the phenotypes)"
)

omxCheckError(
	mxExpectationGREML(V="V",yvars=character(0)),
	"you must specify at least one phenotype in argument 'yvars'"
)

omxCheckError(
	mxExpectationGREML(V="V",yvars="y",yhat="yhat"),
	"non-NULL values for argument 'yhat' are not accepted when 'REML=TRUE'"
)

# In the monophenotype case, we accept a character vector for Xvars:
mxExpectationGREML(V="V",yvars="y",Xvars=c("x1","x2"))
omxCheckError(
	mxExpectationGREML(V="V",yvars=c("y1","y2"),Xvars=c("x1","x2")),
	"argument 'Xvars' must be provided as a list when argument 'yvars' is of length greater than 1"
)

omxCheckError(
	mxExpectationGREML(V="V",yvars="y",Xvars=list(5.5,"foo")),
	"elements of argument 'Xvars' must be of type 'character' (the data column names of the covariates)"
)

omxCheckError(
	mxExpectationGREML(V="V",yvars=c("y1","y2","y3"),Xvars=list("x1","x2",c("x3","x4")),staggerZeroes=FALSE),
	"all phenotypes must have the same number of covariates when staggerZeroes=FALSE"
)

#In the polyphenotype case, the same covariates will often be used for all phenotypes:
mxExpectationGREML(V="V",yvars=c("y1","y2","y3"),Xvars=list(c("x1","x2")))
omxCheckError(
	mxExpectationGREML(V="V",yvars=c("y1","y2","y3"),Xvars=list("x1","x2")),
	"conflicting number of phenotypes specified by arguments 'Xvars' and 'yvars'"
)

omxCheckError(
	mxExpectationGREML(V="V",yvars="y",Xvars="x",REML=FALSE,yhat="yhat"),
	"arguments 'REML' and 'dataset.is.yX' are both FALSE, so non-default values are accepted for only one of arguments 'yhat' (the name of the expected mean vector) or 'Xvars' (the data column names of the covariates, if any)"
)

omxCheckError(
	mxExpectationGREML(V="V",yvars="y",REML=FALSE,yhat=200),
	"argument 'yhat' is not of type 'character' (the name of the expected mean vector)"
)

# In the monophenotype case, we accept a character vector for Xvars:
mxExpectationGREML(V="V",yvars="y",Xvars=c("x1","x2"),REML=FALSE)
omxCheckError(
	mxExpectationGREML(V="V",yvars=c("y1","y2"),Xvars=c("x1","x2"),REML=FALSE),
	"argument 'Xvars' must be provided as a list when argument 'yvars' is of length greater than 1"
)

omxCheckError(
	mxExpectationGREML(V="V",yvars="y",Xvars=list(5.5,"foo"),REML=FALSE),
	"elements of argument 'Xvars' must be of type 'character' (the data column names of the covariates)"
)

omxCheckError(
	mxExpectationGREML(V="V",yvars=c("y1","y2","y3"),Xvars=list("x1","x2",c("x3","x4")),staggerZeroes=FALSE,REML=FALSE),
	"all phenotypes must have the same number of covariates when staggerZeroes=FALSE"
)

#In the polyphenotype case, the same covariates will often be used for all phenotypes:
mxExpectationGREML(V="V",yvars=c("y1","y2","y3"),Xvars=list(c("x1","x2")),REML=FALSE)
omxCheckError(
	mxExpectationGREML(V="V",yvars=c("y1","y2","y3"),Xvars=list("x1","x2"),REML=FALSE),
	"conflicting number of phenotypes specified by arguments 'Xvars' and 'yvars'"
)

testmod <- mxModel(
	"GREMLtest",
	mxData(observed=dat, type="raw", sort=F),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
	mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
	mxAlgebra(I %x% Ve,name="V"),
	mxExpectationGREML(V="V",yvars="y",REML=F,yhat="foo"),
	mxMatrix(type="Full",nrow=100,ncol=1,name="foo",free=T,values=0.12345,labels="bar"),
	mxFitFunctionGREML()
)
testmod@expectation@REML <- TRUE
omxCheckError(
	mxRun(testmod),
	"MxExpecationGREML: slot 'REML' is TRUE and slot 'dataset.is.yX' is FALSE, so slot 'yhat' must have zero length"
)

testmod <- mxModel(
	"GREMLtest",
	mxData(observed=dat, type="raw", sort=F),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
	mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
	mxAlgebra(I %x% Ve,name="V"),
	mxExpectationGREML(V="V",yvars="y",REML=F,yhat="foo"),
	mxMatrix(type="Full",nrow=100,ncol=1,name="foo",free=T,values=0.12345,labels="bar"),
	mxFitFunctionGREML()
)
testmod@expectation@Xvars <- list("x")
omxCheckError(
	mxRun(testmod),
	"MxExpectationGREML: slots 'REML' and 'dataset.is.yX' are both FALSE, so only one of slots 'yhat' and 'Xvars' should have nonzero length"
)

testmod <- mxModel(
	"GREMLtest",
	mxData(observed=dat, type="raw", sort=F),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
	mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
	mxAlgebra(I %x% Ve,name="V"),
	mxExpectationGREML(V="V",yvars="y",REML=F,yhat="foo"),
	mxMatrix(type="Full",nrow=100,ncol=1,name="foo",free=T,values=0.12345,labels="bar"),
	mxMatrix(type="Unit",nrow=100,ncol=1,name="Uno"),
	mxFitFunctionGREML(dyhat=c(ve="Uno",ve="Uno"))
)
omxCheckError(
	mxRun(testmod),
	"duplicated element names in argument 'dyhat'"
)

testmod <- mxModel(
	"GREMLtest",
	mxData(observed=dat, type="raw", sort=F),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
	mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
	mxAlgebra(I %x% Ve,name="V"),
	mxExpectationGREML(V="V",yvars="y",REML=F,yhat="foo"),
	mxMatrix(type="Full",nrow=99,ncol=1,name="foo",free=T,values=0.12345,labels="bar"),
	mxFitFunctionGREML()
)
omxCheckError(
	mxRun(testmod),
	"'y' and 'yhat' vectors have different numbers of rows"
)

testmod <- mxModel(
	"GREMLtest",
	mxData(observed=dat, type="raw", sort=F),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
	mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
	mxAlgebra(I %x% Ve,name="V"),
	mxExpectationGREML(V="V",yvars="y",REML=F,yhat="foo"),
	mxMatrix(type="Full",nrow=100,ncol=1,name="foo",free=T,values=0.12345,labels="bar"),
	mxMatrix(type="Unit",nrow=99,ncol=1,name="Uno"),
	mxFitFunctionGREML(dyhat=c(ve="Uno"))
)
omxCheckError(
	mxRun(testmod),
	"all derivatives of yhat must have the same length as yhat"
)

testmod <- mxModel(
	"GREMLtest",
	mxData(observed=dat, type="raw", sort=F),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
	mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
	mxAlgebra(I %x% Ve,name="V"),
	mxExpectationGREML(V="V",yvars="y",REML=F,yhat="foo"),
	mxMatrix(type="Full",nrow=100,ncol=1,name="foo",free=T,values=0.12345,labels="bar"),
	mxMatrix(type="Unit",nrow=100,ncol=1,name="Uno"),
	mxMatrix(type="Zero",nrow=100,ncol=1,name="Zip"),
	mxAlgebra(Zip %*% t(Zip), name="Zilch"),
	mxFitFunctionGREML(dV=c(ve="I",bar="Zilch",baz="Zilch"),dyhat=c(ve="Zip",bar="Uno",baz="Uno"))
)
omxCheckError(
	mxRun(testmod),
	"length of argument 'dV' is greater than the number of explicit free parameters"
)

testmod$fitfunction <- mxFitFunctionGREML(dV=c(ve="I",bar="Zilch"),dyhat=c(ve="Zip",bar="Uno",baz="Uno"))
omxCheckError(
	mxRun(testmod),
	"length of argument 'dyhat' is greater than the number of explicit free parameters"
)

testmod$fitfunction <- mxFitFunctionGREML(dyhat=c(ve="Zip",bar="Uno",baz="Uno"))
omxCheckError(
	mxRun(testmod),
	"length of argument 'dyhat' is greater than the number of explicit free parameters"
)





# Reset global options: ####
options(mxCondenseMatrixSlots=FALSE)
mxOption(reset=TRUE)
