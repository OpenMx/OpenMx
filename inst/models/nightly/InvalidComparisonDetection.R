
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
#This test does not need to be run with all three GD optimizers:
if(mxOption(NULL,"Default optimizer")!="CSOLNP"){stop("SKIP")}
options(mxCondenseMatrixSlots=TRUE)  #<--Saves memory
require(mvtnorm)


#Generate data:
set.seed(476)
A <- matrix(0,1000,1000)  #<--Empty GRM
A[lower.tri(A)] <- runif(499500, -0.025, 0.025)
A <- A + t(A)
diag(A) <- runif(1000,0.95,1.05) #<--GRM now complete
y <- t(rmvnorm(1,sigma=A*0.5))  #<--Phenotype 'y' has a "population" variance of 1 and h2 of 0.5 
y <- y + rnorm(1000,sd=sqrt(0.5))
x <- rnorm(1000) #<--Covariate 'x' is actually independent of the phenotype.
#Merge variables into data matrix:
dat <- cbind(y,x)
colnames(dat) <- c("y","x") #<--Column names

ge <- mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T)

gff <- mxFitFunctionGREML(dV=c(va="A",ve="I"))

plan <- mxComputeSequence(steps=list(
	mxComputeNewtonRaphson(fitfunction="fitfunction"),
	mxComputeOnce('fitfunction', c('fit','gradient','hessian','ihessian')),
	mxComputeStandardError(),
	mxComputeReportDeriv(),
	mxComputeReportExpectation()
))


mxdat <- mxData(observed = dat, type="raw", sort=FALSE)

testmod <- mxModel(
	"GREML_1GRM_1trait_A", #<--Model name
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = var(y)/2, labels = "ve", lbound = 0.0001, 
					 name = "Ve"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = var(y)/2, labels = "va", name = "Va"),
	mxMatrix("Iden",nrow=1000,name="I"),
	mxMatrix("Symm",nrow=1000,free=F,values=A,name="A"),
	mxAlgebra((A%x%Va) + (I%x%Ve), name="V"),
	mxAlgebra(Va/(Va+Ve), name="h2"),
	mxdat, #<--MxData object
	ge, #<--GREML expectation
	gff, #<--GREML fitfunction
	plan #<--Custom compute plan
)
testrun <- mxRun(testmod)
summary(testrun)

#Drop the covariate:
ge2 <- mxExpectationGREML(V="V",yvars="y", addOnes=T)
testmod2 <- mxModel(
	"GREML_1GRM_1trait_B", #<--Model name
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = var(y)/2, labels = "ve", lbound = 0.0001, 
					 name = "Ve"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = var(y)/2, labels = "va", name = "Va"),
	mxMatrix("Iden",nrow=1000,name="I"),
	mxMatrix("Symm",nrow=1000,free=F,values=A,name="A"),
	mxAlgebra((A%x%Va) + (I%x%Ve), name="V"),
	mxAlgebra(Va/(Va+Ve), name="h2"),
	mxdat, #<--MxData object
	ge2, #<--GREML expectation
	gff, #<--GREML fitfunction
	plan #<--Custom compute plan
)
testrun2 <- mxRun(testmod2)
summary(testrun2)

#Model with fewer free parameters has "better" fit?:
testrun$fitfunction$result
testrun2$fitfunction$result
#But MLfit looks OK:
testrun$fitfunction$MLfit
testrun2$fitfunction$MLfit

omxCheckWarning(
	mxCompare(testrun,testrun2),
	"the names of the fixed effects in MxModels 'GREML_1GRM_1trait_A' and 'GREML_1GRM_1trait_B' do not match; comparison of REML fit values is only valid for models that use the same covariates"
)


rm(testmod,testmod2,testrun,testrun2,A,dat,y,ge,ge2,gff,mxdat,plan,x); gc()


set.seed(1234)
dat <- cbind(rnorm(100),rep(1,100))
colnames(dat) <- c("y","x")

ge <- mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=FALSE)
gff <- mxFitFunctionGREML(dV=c(ve="I"))
plan <- mxComputeSequence(steps=list(
	mxComputeNewtonRaphson(fitfunction="fitfunction"),
	mxComputeOnce('fitfunction', c('fit','gradient','hessian','ihessian')),
	mxComputeStandardError(),
	mxComputeReportDeriv(),
	mxComputeReportExpectation()
))

remlmod <- mxModel(
	"GREMLtest",
	mxData(observed = dat, type="raw", sort=FALSE),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
	mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
	mxAlgebra(I %x% Ve,name="V"),
	ge,
	gff,
	plan
)
remlrun <- mxRun(remlmod)

mlmod <- mxModel(
	"MLtest",
	mxData(observed = dat, type="raw"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0.1,labels="mu",name="Mu"),
	mxExpectationNormal(covariance="Ve",means="Mu",dimnames="y"),
	mxFitFunctionML()
)
mlrun <- mxRun(mlmod)

omxCheckError(
	mxCompare(remlrun,mlrun),
	"MxModel 'GREMLtest' has a fitfunction of class 'MxFitFunctionGREML', but MxModel 'MLtest' has a fitfunction of class 'MxFitFunctionML'"
)


data(factorExample1)

indicators <- names(factorExample1)
latents <- c("F1")
loadingLabels <- paste("b_", indicators, sep="")
uniqueLabels <- paste("U_", indicators, sep="")
meanLabels <- paste("M_", indicators, sep="")
factorVarLabels <- paste("Var_", latents, sep="")

oneFactorCov1 <- mxModel("Single Factor Covariance Model with Fixed Variance",
												 type="RAM",
												 manifestVars=indicators,
												 latentVars=latents,
												 mxPath(from=latents, to=indicators, 
												 			 #           arrows=1, all=TRUE, 
												 			 arrows=1, connect="all.pairs", 
												 			 free=TRUE, values=.2, 
												 			 labels=loadingLabels),
												 mxPath(from=indicators, 
												 			 arrows=2, 
												 			 free=TRUE, values=.8, 
												 			 labels=uniqueLabels),
												 mxPath(from=latents,
												 			 arrows=2, 
												 			 free=FALSE, values=1, 
												 			 labels=factorVarLabels),
												 mxData(observed=cov(factorExample1), type="cov", numObs=500)
)

oneFactorCov1Out <- mxRun(oneFactorCov1)

oneFactorCovWLS <- mxModel(oneFactorCov1Out, name='WLS',
													 mxDataWLS(factorExample1)
)
oneFactorCovWLS <- omxSetParameters(model=oneFactorCovWLS,labels="b_x1",free=F,values=0)

oneFactorCovWLSOut <- mxRun(oneFactorCovWLS)

omxCheckError(
	mxCompare(oneFactorCov1Out,oneFactorCovWLSOut),
	" MxModel 'Single Factor Covariance Model with Fixed Variance' has '-2lnL' fit units, but MxModel 'WLS' has 'r'Wr' fit units"
)
