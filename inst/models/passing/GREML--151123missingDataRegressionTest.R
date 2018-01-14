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
require(mvtnorm)


#Generate data:
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
#y[100] <- NA
x <- rnorm(100) 
dat <- cbind(y,x)
colnames(dat) <- c("y","x")
dat[100,1] <- NA #<--Note that the last row of the dataset is being set to NA.

testmod <- mxModel(
	"GREMLtest",
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values =0.5, labels = "ve", lbound = 0.0001, 
					 name = "Ve"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va1", name = "Va1"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va2", name = "Va2"),
	mxData(observed = dat, type="raw", sort=FALSE),
	mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T),
	mxComputeSequence(steps=list(
		mxComputeNewtonRaphson(fitfunction="fitfunction"),
		mxComputeOnce('fitfunction', c('fit','gradient','hessian','ihessian')),
		mxComputeStandardError(),
		mxComputeReportDeriv()
	)),
	mxMatrix("Iden",nrow=100,name="I"),
	mxMatrix("Symm",nrow=100,free=F,values=A1,name="A1"),
	mxMatrix("Symm",nrow=100,free=F,values=A2,name="A2"),
	mxAlgebra((A1%x%Va1) + (A2%x%Va2) + (I%x%Ve), name="V"),
	mxFitFunctionGREML(dV=c(va1="A1",va2="A2",ve="I"))
)
testrun <- mxRun(testmod) #<--Should run without error.
omxCheckEquals(dim(testrun$V$result),c(100,100))

dat[,1] <- y
dat[90:99,1] <- NA
testmod2 <- mxModel(
	"GREMLtest",
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values =0.5, labels = "ve", lbound = 0.0001, 
					 name = "Ve"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va1", name = "Va1"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va2", name = "Va2"),
	mxData(observed = dat, type="raw", sort=FALSE),
	mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T),
	mxComputeSequence(steps=list(
		mxComputeNewtonRaphson(fitfunction="fitfunction"),
		mxComputeOnce('fitfunction', c('fit','gradient','hessian','ihessian')),
		mxComputeStandardError(),
		mxComputeReportDeriv()
	)),
	mxMatrix("Symm",nrow=100,free=F,values=A1,name="A1"),
	mxMatrix("Symm",nrow=100,free=F,values=A2,name="A2"),
	mxMatrix(type="Unit",nrow=100,ncol=1,name="Uni"),
	#Note that the derivatives of V are all MxAlgebras:
	mxAlgebra(vec2diag(Uni),name="I"),
	mxAlgebra(A1%x%1,name="dV_dva1"),
	mxAlgebra(A2%x%1,name="dV_dva2"),
	mxAlgebra((A1%x%Va1) + (A2%x%Va2) + (I%x%Ve), name="V"),
	mxFitFunctionGREML(dV=c(va1="dV_dva1",va2="dV_dva2",ve="I"))
)
testrun2 <- mxRun(testmod2) #<--Should run without error.
omxCheckEquals(dim(testrun2$I$result),c(100,100))
omxCheckEquals(dim(testrun2$dV_dva1$result),c(100,100))
omxCheckEquals(dim(testrun2$dV_dva2$result),c(100,100))
omxCheckEquals(dim(testrun2$V$result),c(100,100))



#As of commit 384faf7, it should be possible to frontend-downsize derivatives of V that are both MxAlgebras
#and don't depend upon parameters, and things should still work:
testmod3 <- mxModel(
	"GREMLtest",
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values =0.5, labels = "ve", lbound = 0.0001, 
					 name = "Ve"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va1", name = "Va1"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va2", name = "Va2"),
	mxData(observed = dat, type="raw", sort=FALSE),
	mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T),
	mxComputeSequence(steps=list(
		mxComputeNewtonRaphson(fitfunction="fitfunction"),
		mxComputeOnce('fitfunction', c('fit','gradient','hessian','ihessian')),
		mxComputeStandardError(),
		mxComputeReportDeriv()
	)),
	#The two GRMs right below are dependencies of V, so they need to be 100x100:
	mxMatrix("Symm",nrow=100,free=F,values=A1,name="A1"),
	mxMatrix("Symm",nrow=100,free=F,values=A2,name="A2"),
	mxMatrix(type="Unit",nrow=100,ncol=1,name="Uni"),
	#Derivatives of V:
	mxMatrix("Symm",nrow=90,free=F,values=A1[c(1:89,100),c(1:89,100)],name="dV_dva1"),
	mxAlgebra(A2%x%1,name="dV_dva2"),
	mxAlgebra(vec2diag(Uni[1:90,1]), name="dV_dve"),
	mxAlgebra((A1%x%Va1) + (A2%x%Va2) + vec2diag(Uni%x%Ve), name="V"),
	mxFitFunctionGREML(dV=c(va1="dV_dva1",va2="dV_dva2",ve="dV_dve"))
)
testrun3 <- mxRun(testmod3)
#The important thing is that everything come back with the dimensions the user has said they should be:
omxCheckEquals(dim(testrun3$V$result),c(100,100))
omxCheckEquals(dim(testrun3$dV_dva1$values),c(90,90))
omxCheckEquals(dim(testrun3$dV_dva2$result),c(100,100))
omxCheckEquals(dim(testrun3$dV_dve$result),c(90,90))



#For a simple model like this, it should also be possible to frontend-downsize everything, if you call
#the data-handler before runtime and plan accordingly:
gremldat <- mxGREMLDataHandler(data=dat,yvars="y",Xvars="x")
testmod4 <- mxModel(
	"GREMLtest",
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values =0.5, labels = "ve", lbound = 0.0001, 
					 name = "Ve"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va1", name = "Va1"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va2", name = "Va2"),
	mxData(observed = gremldat$yX, type="raw", sort=FALSE),
	mxExpectationGREML(V="V",dataset.is.yX=TRUE), #<--Note we don't provide any casesToDropFromV.
	mxComputeSequence(steps=list(
		mxComputeNewtonRaphson(fitfunction="fitfunction"),
		mxComputeOnce('fitfunction', c('fit','gradient','hessian','ihessian')),
		mxComputeStandardError(),
		mxComputeReportDeriv()
	)),
	mxMatrix("Symm",nrow=90,free=F,values=A1[c(1:89,100),c(1:89,100)],name="A1"),
	mxMatrix("Symm",nrow=90,free=F,values=A2[c(1:89,100),c(1:89,100)],name="A2"),
	mxMatrix("Iden",nrow=90,name="I"),
	mxAlgebra((A1%x%Va1) + (A2%x%Va2) + I%x%Ve, name="V"),
	mxFitFunctionGREML(dV=c(va1="A1",va2="A2",ve="I"))
)
testrun4 <- mxRun(testmod4) #<--Should run without error.


