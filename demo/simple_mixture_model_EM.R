#
#   Copyright 2007-2024 by the individuals mentioned in the source code history
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

# -----------------------------------------------------------------------
# Script by Robert M. Kirkpatrick.
# OpenMx Master Class, November 2024, Richmond, Virginia, USA.

# This script demonstrates fitting a univariate, 2-class mixture model
# via expectation-maximization (EM) in OpenMx.

# Load OpenMx: ####
library(OpenMx)

# Simulate data: ####
set.seed(190127) #<--Set RNG seed:
mu1 <- -1.5 #<--Class 1's true mean.
mu2 <- 1.5 #<--Class 2's true mean.
N <- 200 #<--Total sample size.
# Half the sample is from class 1, and the other half is from class 2;
# both classes happen to have the same variance:
x <- matrix(
	c(rnorm(N/2,mu1,1),
	rnorm(N/2,mu2,1)),ncol=1,dimnames=list(NULL,"x"))
# MxData object:
data4mx <- mxData(observed=x,type="raw")

# MxModel object for class 1: ####
class1 <- mxModel(
	"Class1", #<--Model name.
	mxMatrix(type="Full",nrow=1,ncol=1,free=TRUE,values=-.5,name="Mu"), #<--Mean.
	mxMatrix(type="Full",nrow=1,ncol=1,free=TRUE,values=4,name="Sigma"), #<--Variance.
	#^^^N.B. the two MxMatrices above have been given names, but their elements have
	# not been given parameter labels!  In the OpenMx namespace, labels are global,
	# but names are local.
	mxExpectationNormal(covariance="Sigma",means="Mu",dimnames="x"), #<--Expectation.
	# Because of `vector=TRUE`, the fitfunction will return a vector of raw likelihoods
	# (it usually returns a scalar sum of -2loglikelihoods):
	mxFitFunctionML(vector=TRUE))

# MxModel object for class 2: ####
class2 <- mxRename(class1, "Class2")
class2$Mu$values[1,1] <- .5

# Mixture model: ####
mm <- mxModel(
	"Mixture",  #<--Model name.
	data4mx, #<--MxData object.
	class1, class2, #<--Child MxModels for the two mixture classes.
	# Class-1 likelihoods weighted by probability of class-1 membership:
	mxAlgebra((1-Posteriors) * Class1.fitfunction, name="PL1"),
	# Class-2 likelihoods weighted by probability of class-2 membership:
	mxAlgebra(Posteriors * Class2.fitfunction, name="PL2"),
	# Sums of the 2 classes' probability-weighted likelihoods:
	mxAlgebra(PL1 + PL2, name="PL"),
	# Posterior probabilities of class-2 membership:
	mxAlgebra(
		PL2 / PL,  recompute='onDemand',
		initial=matrix(runif(N,.4,.6), nrow=N, ncol = 1), name="Posteriors"),
	#^^^N.B. `recompute` & `initial`!
	mxAlgebra(-2*sum(log(PL)), name="FF"), #<-- -2logL for whole sample.
	mxFitFunctionAlgebra(algebra="FF"), #<--Algebra fitfunction.
	# The custom compute plan:
	mxComputeEM(
		estep=mxComputeOnce("Mixture.Posteriors"),
		mstep=mxComputeGradientDescent(fitfunction="Mixture.fitfunction")))

# Run the model & see output: ####
mmfit <- mxRun(mm)
( smm <- summary(mmfit) )
mxEval(Posteriors,mmfit) #<--Posterior probabilities of class-2 membership.

omxCheckCloseEnough(coef(mmfit),c(-1.306432,1.024295,1.884625,0.716655),0.001)
omxCheckCloseEnough(mmfit$output$fit,539.1599,0.0001)
omxCheckEquals(smm$degreesOfFreedom,196)
