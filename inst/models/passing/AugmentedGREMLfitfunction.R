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

#Baseline model:
testmod <- mxModel(
	"GREMLtest",
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values =0.5, labels = "ve", lbound = 0.0001, 
					 name = "Ve"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va1", name = "Va1"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va2", name = "Va2"),
	mxData(observed = dat, type="raw", sort=FALSE),
	mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T),
	mxMatrix("Iden",nrow=100,name="I"),
	mxMatrix("Symm",nrow=100,free=F,values=A1,name="A1"),
	mxMatrix("Symm",nrow=100,free=F,values=A2,name="A2"),
	mxAlgebra((A1%x%Va1) + (A2%x%Va2) + (I%x%Ve), name="V"),
	mxFitFunctionGREML()
)
testrun <- mxRun(testmod)

#Pointless augmentation that adds a constant to the fitfunction:
testmod2 <- mxModel(
	"GREMLtest",
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values =0.5, labels = "ve", lbound = 0.0001, 
					 name = "Ve"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va1", name = "Va1"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va2", name = "Va2"),
	mxData(observed = dat, type="raw", sort=FALSE),
	mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T),
	mxMatrix("Iden",nrow=100,name="I"),
	mxMatrix("Symm",nrow=100,free=F,values=A1,name="A1"),
	mxMatrix("Symm",nrow=100,free=F,values=A2,name="A2"),
	mxAlgebra((A1%x%Va1) + (A2%x%Va2) + (I%x%Ve), name="V"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=F,values=0.64,name="aug"),
	mxFitFunctionGREML(aug="aug")
)
testrun2 <- mxRun(testmod2)
omxCheckCloseEnough(a=testrun2$output$fit - testrun$output$fit, b=1.28, epsilon=1e-9)

#Baseline model using N-R:
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
		mxComputeReportDeriv(),
		mxComputeReportExpectation()
	)),
	mxMatrix("Iden",nrow=100,name="I"),
	mxMatrix("Symm",nrow=100,free=F,values=A1,name="A1"),
	mxMatrix("Symm",nrow=100,free=F,values=A2,name="A2"),
	mxAlgebra((A1%x%Va1) + (A2%x%Va2) + (I%x%Ve), name="V"),
	mxFitFunctionGREML(dV=c(va1="A1",va2="A2",ve="I"))
)
testrun3 <- mxRun(testmod3)

#Add augmentation that should nudge free parameters toward summing to 1.0:
testmod4 <- mxModel(
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
		mxComputeReportDeriv(),
		mxComputeReportExpectation()
	)),
	mxMatrix("Iden",nrow=100,name="I"),
	mxMatrix("Symm",nrow=100,free=F,values=A1,name="A1"),
	mxMatrix("Symm",nrow=100,free=F,values=A2,name="A2"),
	mxAlgebra((A1%x%Va1) + (A2%x%Va2) + (I%x%Ve), name="V"),
	mxAlgebra( 3%x%(Va1+Va2+Ve-1)^2, name="aug"),
	mxAlgebra( 3%x%rbind(
		2*Va1 + 2*Va2 + 2*Ve - 2,
		2*Va1 + 2*Va2 + 2*Ve - 2,
		2*Va1 + 2*Va2 + 2*Ve - 2), name="daug1"),
	mxMatrix(type="Full",nrow=3,ncol=3,free=F,values=6,name="daug2"),
	mxFitFunctionGREML(dV=c(va1="A1",va2="A2",ve="I"),aug="aug",augGrad="daug1",augHess="daug2")
)
testrun4 <- mxRun(testmod4)
#The difference between 1.0 and the sum of the parameters should be smaller for model #4:
omxCheckTrue(abs(1-sum(testrun4$output$estimate)) < abs(1-sum(testrun3$output$estimate)))

