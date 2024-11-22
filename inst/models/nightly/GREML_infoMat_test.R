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

require(OpenMx)
#This test doesn't do any optimization, so there's no need to run it with all 3
#gradient-based optimizers:
if(mxOption(NULL,"Default optimizer")!="SLSQP"){stop("SKIP")}
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

#The GREML expectation tells OpenMx that the model-expected covariance matrix is named 'V', that the one 
#phenotype is has column label 'y' in the dataset, that the one covariate has column label 'x' in the dataset,
#and that a lead column of ones needs to be appended to the 'X' matrix (for the intercept):
ge <- mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T)

#This is a custom compute plan.  It is necessary because we want to use the Newton-Raphson optimizer, which
#can use analytic first and second derivatives of the GREML fitfunction to speed up convergence.  It looks
#especially messy here because we want a profile-likelihood confidence interval for the heritability:
plan <- mxComputeSequence(steps=list(
	mxComputeOnce('fitfunction', c('fit','gradient','hessian')),
	mxComputeReportDeriv(),
	mxComputeReportExpectation()
))

#The MxData object.  N.B. use of 'sort=FALSE' is CRITICALLY IMPORTANT, because the rows and columns of dataset
#'dat' and the rows and columns of GRM 'A' are already properly ordered:
mxdat <- mxData(observed = dat, type="raw", sort=FALSE)

#We will create some of the necessary objects inside the mxModel() statement.  We mainly want to avoid creating 
#more copies of the GRM than we need to:
testmod <- mxModel(
	"GREML_1GRM_1trait", #<--Model name
	#1x1 matrix containing residual variance component:
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = var(y)/2, labels = "ve", lbound = 0.0001, 
					 name = "Ve"),
	#1x1 matrix containing additive-genetic variance component:
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = var(y)/2, labels = "va", name = "Va"),
	#1000x1000 identity matrix--the "relatedness matrix" for the residuals:
	mxMatrix("Iden",nrow=1000,name="I"),
	#The GRM:
	mxMatrix("Symm",nrow=1000,free=F,values=A,name="A"),
	#The model-expected covariance matrix:
	mxAlgebra((A%x%Va) + (I%x%Ve), name="V"),
	#An MxAlgebra for the heritability:
	mxAlgebra(Va/(Va+Ve), name="h2"),
	mxCI("h2"), #<--Request confidence interval for heritability
	mxdat, #<--MxData object
	ge, #<--GREML expectation
	mxFitFunctionGREML(dV=c(va="A",ve="I"),infoMatType="expected"), #<--GREML fitfunction
	plan #<--Custom compute plan
)

testrun <- mxRun(testmod) #<--Run model
testrun$output$gradient
testrun$output$hessian

mygrad <- c(ve=NA_real_,va=NA_real_)
myEIM <- matrix(NA_real_,2,2,dimnames=list(c("ve","va"),c("ve","va")))

Vinv <- mxEval(chol2inv(chol(V)),testrun,T)
dat2 <- mxGREMLDataHandler(data=dat,yvars="y", Xvars="x", addOnes=T)
X <- dat2$yX[,-1]
y <- as.matrix(dat2$yX[,1])
P <- Vinv - Vinv%*%X%*%chol2inv(chol(t(X)%*%Vinv%*%X))%*%t(X)%*%Vinv
I <- diag(1000)
mygrad[1] <- tr(I%*%P) - t(P%*%y)%*%I%*%P%*%y
mygrad[2] <- tr(A%*%P) - t(P%*%y)%*%A%*%P%*%y
omxCheckCloseEnough(mygrad, testrun$output$gradient, 5e-12)

myEIM[1,1] <- sum((I%*%P)*t(I%*%P))
myEIM[1,2] <- sum((I%*%P)*t(A%*%P) )
myEIM[2,1] <- myEIM[1,2]
myEIM[2,2] <- sum((A%*%P)*t(A%*%P))
omxCheckCloseEnough(myEIM, testrun$output$hessian, 1e-10)


testmod2 <- mxModel(
	"GREML_1GRM_1trait", #<--Model name
	#1x1 matrix containing residual variance component:
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = var(y)/2, labels = "ve", lbound = 0.0001, 
					 name = "Ve"),
	#1x1 matrix containing additive-genetic variance component:
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = var(y)/2, labels = "va", name = "Va"),
	#1000x1000 identity matrix--the "relatedness matrix" for the residuals:
	mxMatrix("Iden",nrow=1000,name="I"),
	#The GRM:
	mxMatrix("Symm",nrow=1000,free=F,values=A,name="A"),
	#The model-expected covariance matrix:
	mxAlgebra((A%x%Va) + (I%x%Ve), name="V"),
	#An MxAlgebra for the heritability:
	mxAlgebra(Va/(Va+Ve), name="h2"),
	mxCI("h2"), #<--Request confidence interval for heritability
	mxdat, #<--MxData object
	ge, #<--GREML expectation
	mxFitFunctionGREML(dV=c(va="A",ve="I"),infoMatType="average"), #<--GREML fitfunction
	plan #<--Custom compute plan
)

testrun2 <- mxRun(testmod2) #<--Run model
testrun2$output$gradient
testrun2$output$hessian

mygrad[1] <- sum(I*P) - t(P%*%y)%*%I%*%P%*%y
mygrad[2] <- sum(A*P) - t(P%*%y)%*%A%*%P%*%y
omxCheckCloseEnough(mygrad, testrun2$output$gradient, 5e-12)

myEIM[1,1] <- t(P%*%y)%*%I%*%P%*%I%*%P%*%y
myEIM[1,2] <- t(P%*%y)%*%I%*%P%*%A%*%P%*%y
myEIM[2,1] <- myEIM[1,2]
myEIM[2,2] <- t(P%*%y)%*%A%*%P%*%A%*%P%*%y
omxCheckCloseEnough(myEIM, testrun2$output$hessian, 1e-10)

options(mxCondenseMatrixSlots=FALSE)
