#
#   Copyright 2007-2025 by the individuals mentioned in the source code history
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
# mxOption(NULL,"Default optimizer","NPSOL")
# mxOption(NULL,"Number of threads",2)
# mxOption(NULL,"Print level",20)
# mxOption(NULL,"Print file",3)
# mxOption(NULL,"Verify level",3)
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
y[100] <- NA
x <- rnorm(100)
dat <- cbind(y,x)
colnames(dat) <- c("y","x")

plan <- mxComputeSequence(
	steps=list(
		mxComputeGradientDescent(engine=mxOption(NULL,"Default optimizer")),
		mxComputeOnce('fitfunction', c('gradient','hessian')),
		mxComputeStandardError(),
		mxComputeHessianQuality(),
		mxComputeReportDeriv(),
		mxComputeReportExpectation()
	))

for(pds in 1:3){
	test0 <- mxModel(
		"GREMLtest0",
		plan,
		mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values =0.5, labels = "ve", lbound = 0.0001,
						 name = "Ve"),
		mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va1", name = "Va1"),
		mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.20, labels = "va2", name = "Va2"),
		mxData(observed = dat, type="raw", sort=FALSE),
		mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T),
		mxMatrix("Iden",nrow=100,name="I"),
		mxMatrix("Symm",nrow=100,free=F,values=A1,name="A1"),
		mxMatrix("Symm",nrow=100,free=F,values=A2,name="A2"),
		mxAlgebra((A1%x%Va1) + (A2%x%Va2) + (I%x%Ve), name="V"),
		mxFitFunctionGREML(dV=c(va1="A1",va2="A2",ve="I"))
	)
	test0$fitfunction@.parallelDerivScheme <- pds
	test0 <- mxRun(test0)
	
	test1 <- mxModel(
		"GREMLtest1",
		plan,
		mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values =0.5, labels = "ve", lbound = 0.0001,
						 name = "Ve"),
		mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va1", name = "Va1"),
		mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.20, labels = "va2", name = "Va2"),
		mxData(observed = dat, type="raw", sort=FALSE),
		mxExpectationGREML(V="V",yvars="y", Xvars=list(), addOnes=T),
		mxMatrix("Iden",nrow=100,name="I"),
		mxMatrix("Symm",nrow=100,free=F,values=A1,name="A1"),
		mxMatrix("Symm",nrow=100,free=F,values=A2,name="A2"),
		mxAlgebra((A1%x%Va1) + (A2%x%Va2) + (I%x%Ve), name="V"),
		mxFitFunctionGREML(dV=c(va1="A1",va2="A2",ve="I"))
	)
	test1$fitfunction@.parallelDerivScheme <- pds
	test1 <- mxRun(test1)
	
	test2 <- mxModel(
		"GREMLtest2",
		plan,
		mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values =0.5, labels = "ve", lbound = 0.0001,
						 name = "Ve"),
		mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va1", name = "Va1"),
		mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.20, labels = "va2", name = "Va2"),
		mxData(observed = dat, type="raw", sort=FALSE),
		mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T, REML=F),
		mxMatrix("Iden",nrow=100,name="I"),
		mxMatrix("Symm",nrow=100,free=F,values=A1,name="A1"),
		mxMatrix("Symm",nrow=100,free=F,values=A2,name="A2"),
		mxAlgebra((A1%x%Va1) + (A2%x%Va2) + (I%x%Ve), name="V"),
		mxFitFunctionGREML(dV=c(va1="A1",va2="A2",ve="I"))
	)
	test2$fitfunction@.parallelDerivScheme <- pds
	test2 <- mxRun(test2)
	
	test3 <- mxModel(
		"GREMLtest3",
		plan,
		mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values =0.5, labels = "ve", lbound = 0.0001,
						 name = "Ve"),
		mxMatrix(type = "Full", nrow = 1, ncol=1, free=F, values = 0, labels = "va1", name = "Va1"),
		mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.20, labels = "va2", name = "Va2"),
		mxData(observed = dat, type="raw", sort=FALSE),
		mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T),
		mxMatrix("Iden",nrow=100,name="I"),
		mxMatrix("Symm",nrow=100,free=F,values=A2,name="A2"),
		mxAlgebra((A2%x%Va2) + (I%x%Ve), name="V"),
		mxFitFunctionGREML(dV=c(va2="A2",ve="I"))
	)
	test3$fitfunction@.parallelDerivScheme <- pds
	test3 <- mxRun(test3)
	
	test4 <- mxModel(
		"GREMLtest4",
		plan,
		mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values =0.5, labels = "ve", lbound = 0.0001,
						 name = "Ve"),
		mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va1", name = "Va1"),
		mxMatrix(type = "Full", nrow = 1, ncol=1, free=F, values = 0, labels = "va2", name = "Va2"),
		mxData(observed = dat, type="raw", sort=FALSE),
		mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T),
		mxMatrix("Iden",nrow=100,name="I"),
		mxMatrix("Symm",nrow=100,free=F,values=A1,name="A1"),
		mxAlgebra((A1%x%Va1) + (I%x%Ve), name="V"),
		mxFitFunctionGREML(dV=c(va1="A1",ve="I"))
	)
	test4$fitfunction@.parallelDerivScheme <- pds
	test4 <- mxRun(test4)
	
	test5 <- mxModel(
		"GREMLtest5",
		plan,
		mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values =0.5, labels = "ve", lbound = 0.0001,
						 name = "Ve"),
		mxMatrix(type = "Full", nrow = 1, ncol=1, free=F, values = 0, labels = "va1", name = "Va1"),
		mxMatrix(type = "Full", nrow = 1, ncol=1, free=F, values = 0, labels = "va2", name = "Va2"),
		mxData(observed = dat, type="raw", sort=FALSE),
		mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T),
		mxMatrix("Iden",nrow=100,name="I"),
		mxAlgebra(I%x%Ve, name="V"),
		mxFitFunctionGREML(dV=c(ve="I"))
	)
	test5$fitfunction@.parallelDerivScheme <- pds
	test5 <- mxRun(test5)
	
	mxModelAverage(reference=c("ve","va1","va2"),models=list(test0,test3,test4,test5)) #<--Should work.
	omxCheckWarning(
		mxModelAverage(reference=c("ve","va1","va2"),models=list(test0,test1,test3,test4,test5)),
		"not all of 'models' have matching names for their fixed effects; results may be invalid"
	)
	omxCheckError(
		mxModelAverage(reference=c("ve","va1","va2"),models=list(test0,test1,test2,test3,test4,test5)),
		"some but not all of 'models' use REML"
	)
	
	test6 <- mxModel(
		"GREMLtest6",
		plan,
		mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values =0.5, labels = "ve", lbound = 0.0001,
						 name = "Ve"),
		mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va1", name = "Va1"),
		mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.20, labels = "va2", name = "Va2"),
		mxData(observed = dat, type="raw", sort=FALSE),
		mxExpectationGREML(V="V",yvars="y", Xvars=list(), addOnes=T, REML=F),
		mxMatrix("Iden",nrow=100,name="I"),
		mxMatrix("Symm",nrow=100,free=F,values=A1,name="A1"),
		mxMatrix("Symm",nrow=100,free=F,values=A2,name="A2"),
		mxAlgebra((A1%x%Va1) + (A2%x%Va2) + (I%x%Ve), name="V"),
		mxFitFunctionGREML(dV=c(va1="A1",va2="A2",ve="I"))
	)
	test6$fitfunction@.parallelDerivScheme <- pds
	test6 <- mxRun(test6)
	mxModelAverage(reference=c("ve","va1","va2"),models=list(test2,test6)) #<--Should work.
}




options(mxCondenseMatrixSlots=FALSE)
