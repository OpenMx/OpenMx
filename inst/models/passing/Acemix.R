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


# -----------------------------------------------------------------------
# Mixture distribution model: probabilistically diagnosed zygosity
#
# After Neale (2003) A finite mixture distribution model for data collected from twins. 
#       Twin Res 2003; 6:235-239
# ACE Model is specified with RawData and Matrix-style Input
#
# Uses two probabilities: that of correct zygosity diagnosis, pright; 
# and that of incorrect, pwrong. These could be different for MZ and DZ pairs
# but are not here. Using definition variables each pair could have its own
# probability of correct diagnosis. That is not done here either.
# -----------------------------------------------------------------------

require(OpenMx)
dataMZ <-suppressWarnings(try(read.table("data/sim1.mz",header=F), silent=TRUE))
if (is(dataMZ, "try-error")) dataMZ <- read.table("models/passing/data/sim1.mz",header=F)
dataDZ <-suppressWarnings(try(read.table("data/sim1.dz",header=F), silent=TRUE))
if (is(dataDZ, "try-error")) dataDZ <- read.table("models/passing/data/sim1.dz",header=F)
dataMZ <- dataMZ[,1:2] # discard unused data column
dataDZ <- dataDZ[,1:2]
selVars<-c("T1","T2")
names(dataMZ)<-selVars
names(dataDZ)<-selVars
pright <- .95
pwrong <- 1 - pright

twinACEModel <- mxModel("twinACE", 
	# Matrix expMean for expected mean vector for MZ and DZ twins    
	mxMatrix("Full", 1, 2, T, 20, "mean", name="expMean"), 
	# Matrices X, Y, and Z to store the a, c, and e path coefficients
	mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=.6, lbound=.01 , label="a", name="X"),
	mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=.6, lbound=.01 , label="c", name="Y"),
	mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=.6, lbound=.1, ubound=10, label="e", name="Z"),
	# Matrixes A, C, and E to compute A, C, and E variance components
	mxAlgebra(X * t(X), name="A"),
	mxAlgebra(Y * t(Y), name="C"),
	mxAlgebra(Z * t(Z), name="E"),  

	# Matrix expCOVMZ for expected covariance matrix for MZ twins
	mxAlgebra(rbind(cbind(A+C+E   , A+C),
	                cbind(A+C     , A+C+E)), name="expCovMZ"),

	# Matrix expCOVMZ for expected covariance matrix for DZ twins
	mxAlgebra(rbind(cbind(A+C+E   , .5%x%A+C),
	                cbind(.5%x%A+C , A+C+E)), name="expCovDZ"),
	# MZ likelihood is set up as pright*(Likelihood|Zygosity=MZ) + pwrong*(Likelihood|DZ)
	# DZ likelihood is set up as pright*(Likelihood|Zygosity=DZ) + pwrong*(Likelihood|MZ)
	# vector=TRUE allows mixture distribution on individual likelihoods
    mxModel("MZcorrect",
            mxData(dataMZ, type="raw"),
            mxExpectationNormal("twinACE.expCovMZ", "twinACE.expMean",selVars),
            mxFitFunctionML(vector=T)),
    mxModel("MZincorrect",
            mxData(dataMZ, type="raw"),
            mxExpectationNormal("twinACE.expCovDZ", "twinACE.expMean",selVars),
            mxFitFunctionML(vector=T)),
    mxModel("DZcorrect", 
            mxData(dataDZ, type="raw"), 
            mxExpectationNormal("twinACE.expCovDZ", "twinACE.expMean",selVars),
            mxFitFunctionML(vector=T)),
    mxModel("DZincorrect", 
            mxData(dataDZ, type="raw"), 
            mxExpectationNormal("twinACE.expCovMZ", "twinACE.expMean",selVars),
            mxFitFunctionML(vector=T)),
    mxAlgebra(-2*sum(log(pright%x%MZcorrect.fitfunction + pwrong%x%MZincorrect.fitfunction)) + 
              -2*sum(log(pright%x%DZcorrect.fitfunction + pwrong%x%DZincorrect.fitfunction)), name="twin"), 
    mxFitFunctionAlgebra("twin")
)
twinACEFit <- mxRun(twinACEModel, suppressWarnings=TRUE)

# Check results against hard-coded Mx1 estimates and likelihood
estimates     <- mxEval(as.vector(c(X, Y, Z, expMean[1,1])), twinACEFit)
fitStatistics <- mxEval(fitfunction, twinACEFit)
omxCheckCloseEnough(as.vector(c(0.7974,0.4196,0.4509,-0.0254,10165.966)),as.vector(c(estimates,fitStatistics)), .005)
