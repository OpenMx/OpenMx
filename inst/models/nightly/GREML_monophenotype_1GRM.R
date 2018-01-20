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

#This script (by Rob K.) demonstrates the use of the GREML feature in a simple but realistic example.
#It first simulates a genomic-relatedness matrix (GRM), a phenotype, and a null covariate.  Then, it
#fits a simple GREML model to estimate additive-genetic variance, residual variance, and heritability.


require(OpenMx)
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

#The GREML fitfunction tells OpenMx that the derivative of 'V' with respect to free parameter 
#'va'(the additive-genetic variance) is a matrix named 'A', and that the derivative of 'V' w/r/t free parameter
#'#'ve' is a matrix named 'I'.  At runtime, the GREML fitfunction will use these derivatives to help with 
#'#optimization and compute standard errors:
gff <- mxFitFunctionGREML(dV=c(va="A",ve="I"))

#This is a custom compute plan.  It is necessary because we want to use the Newton-Raphson optimizer, which
#can use analytic first and second derivatives of the GREML fitfunction to speed up convergence.  It looks
#especially messy here because we want a profile-likelihood confidence interval for the heritability:
plan <- mxComputeSequence(steps=list(
	mxComputeNewtonRaphson(fitfunction="fitfunction"),
	mxComputeOnce('fitfunction', c('fit','gradient','hessian','ihessian')),
	mxComputeConfidenceInterval(
		plan=mxComputeGradientDescent(
			fitfunction="GREML_1GRM_1trait.fitfunction", nudgeZeroStarts=FALSE, maxMajorIter=150),
		fitfunction="GREML_1GRM_1trait.fitfunction"),
	mxComputeStandardError(),
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
	gff, #<--GREML fitfunction
	plan #<--Custom compute plan
)

testrun <- mxRun(testmod,intervals = T) #<--Run model (Status Red is OK in this case)
summary(testrun) #<--Model summary

#Obtain SE of h2 from delta-method approximation (e.g., Lynch & Walsh, 1998, Appendix 1):
scm <- chol2inv(chol(testrun$output$hessian/2)) #<--Sampling covariance matrix for ve and va
pointest <- testrun$output$estimate #<--Point estimates of ve and va
h2se <- sqrt(
	(pointest[2]/(pointest[1]+pointest[2]))^2 * (
		(scm[2,2]/pointest[2]^2) - (2*scm[1,2]/pointest[1]/(pointest[1]+pointest[2])) + 
			(sum(scm)*(pointest[1]+pointest[2])^-2)
))
#Compare:
mxEval(h2,testrun,T) + 2*c(-h2se,h2se)
testrun$output$confidenceIntervals

#Test for regressions in how GREML handles its analytic derivatives (for OpenMx developer use):
omxCheckTrue(testrun$output$hessian[1,2]>0)
omxCheckCloseEnough(testrun$output$gradient,c(0,0),epsilon=0.15)



#Diagonalize the problem: ###
eigenA <- eigen(A) #<--Eigen decomposition of the GRM
#We "rotate out" the dependence among participants by premultiplying the 'y' and 'X' matrices by the 
#eigenvectors of the GRM:
yrot <- t(eigenA$vectors) %*% y
xrot <- t(eigenA$vectors) %*% cbind(1,x)
datrot <- cbind(yrot,xrot)
colnames(datrot) <- c("y","x0","x1")
#Make a new MxModel:
testmod2 <- mxModel(
	"GREMLtest_1GRM_1trait_diagonalized",
	mxData(observed = datrot, type="raw", sort=FALSE),
	mxExpectationGREML(V="V",yvars="y", Xvars=c("x0","x1"), addOnes=F),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = var(y)/2, labels = "ve", lbound = 0.0001, 
					 name = "Ve"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = var(y)/2, labels = "va", name = "Va"),
	mxMatrix("Iden",nrow=1000,name="I"),
	mxMatrix("Diag",nrow=1000,free=F,values=eigenA$values,name="A"),
	mxAlgebra((A%x%Va) + (I%x%Ve), name="V"),
	mxAlgebra(Va/(Va+Ve), name="h2"),
	gff,
	#We'll do without the CI this time:
	mxComputeSequence(steps=list(
		mxComputeNewtonRaphson(fitfunction="fitfunction"),
		mxComputeOnce('fitfunction', c('fit','gradient','hessian','ihessian')),
		mxComputeStandardError(),
		mxComputeReportDeriv(),
		mxComputeReportExpectation()
	))
)

testrun2 <- mxRun(testmod2)
#Results are substantially equivalent to those from the previous MxModel:
summary(testrun2)
mxEval(h2,testrun2,T)

#OpenMx developer tests:
omxCheckCloseEnough(testrun2$output$gradient,c(0,0),1e-3)
omxCheckCloseEnough(testrun$output$hessian,testrun2$output$hessian,epsilon=1)
omxCheckCloseEnough(testrun$output$estimate,testrun2$output$estimate,0.001)
omxCheckCloseEnough(testrun$output$standardErrors,testrun2$output$standardErrors,0.0001)
omxCheckCloseEnough(testrun$expectation$b,testrun2$expectation$b,1e-5)
omxCheckCloseEnough(testrun$expectation$bcov,testrun2$expectation$bcov,1e-6)


# Reparametrize the problem in terms of total variance and heritability: ###

gff3 <- mxFitFunctionGREML(dV=c(h2="dVdH2",vp="dVdVp")) #<--Need new fitfunction object
#Need new compute plan:
plan3 <- mxComputeSequence(steps=list(
	mxComputeNewtonRaphson(fitfunction="fitfunction"),
	mxComputeOnce('fitfunction', c('fit','gradient','hessian','ihessian')),
	mxComputeConfidenceInterval(
		plan=mxComputeGradientDescent(
			fitfunction="GREML_1GRM_1trait_altparam.fitfunction", nudgeZeroStarts=FALSE, maxMajorIter=150),
		fitfunction="GREML_1GRM_1trait_altparam.fitfunction"),
	mxComputeStandardError(),
	mxComputeReportDeriv(),
	mxComputeReportExpectation()
))

testmod3 <- mxModel(
	"GREML_1GRM_1trait_altparam", #<--Model name
	#1x1 matrix containing heritability:
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.5, labels = "h2", lbound = 0.0001, ubound=0.9999,
					 name = "H2"),
	#1x1 matrix containing total phenotypic variance:
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = var(y), labels = "vp", name = "Vp"),
	#1000x1000 identity matrix--the "relatedness matrix" for the residuals:
	mxMatrix("Iden",nrow=1000,name="I"),
	#The GRM:
	mxMatrix("Symm",nrow=1000,free=F,values=A,name="A"),
	#MxAlgebra for additive-genetic variance:
	mxAlgebra(H2*Vp, name="Va"),
	#MxAlgebra for residual variance:
	mxAlgebra((1-H2)*Vp, name="Ve"),
	#The model-expected covariance matrix:
	mxAlgebra( (A%x%Va) + (I%x%Ve), name="V"),
	#MxAlgebras for derivatives of V w/r/t free parameters:
	mxAlgebra((A-I)%x%Vp, name="dVdH2"),
	mxAlgebra(I + (A-I)%x%H2, name="dVdVp"),
	mxCI("h2"), #<--Request confidence interval for heritability
	mxdat, #<--MxData object
	ge, #<--GREML expectation
	gff3, #<--GREML fitfunction
	plan3 #<--Custom compute plan
)

testrun3 <- mxRun(testmod3, intervals = T)
summary(testrun3)

#Compare:
mxEval(h2,testrun3,T) + 2*c(-0.07824315,0.07824315) #<--0.07824315 is the SE of h2
testrun3$output$confidenceIntervals
