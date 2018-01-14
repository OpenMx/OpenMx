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
#
# ACE Model is specified with RawData and Matrix-style Input
#
# Uses two probabilities: that of correct zygosity diagnosis, pright; 
# and that of incorrect, pwrong.  These could be different for MZ and DZ pairs
# but are not here.  Using definition variables each pair could have its own
# probability of correct diagnosis.  That is not done here either.
#
# -----------------------------------------------------------------------

require(OpenMx)

DataMZ <- suppressWarnings(try(read.table("models/passing/data/sim1.mz", header=FALSE), silent=TRUE))
if (is(DataMZ, "try-error")) DataMZ <- read.table("data/sim1.mz", header = F)
selVars <- c("T1", "T2", "pMZ")
names(DataMZ) <- selVars
frameMZ <- data.frame(pMZ = DataMZ$pMZ)

twinACEModel <- mxModel("twinACE", 
	        # Matrix expMean for expected mean vector for MZ and DZ twins    
	mxMatrix("Full", nrow=1, ncol=2, free=TRUE, values=0,
	         label="mean", name="expMean"), 

	        # Matrices X, Y, and Z to store the a, c, and e path coefficients
	mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=.6, label="a", name="X"),
	mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=.6, label="c", name="Y"),
	mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=.6, label="e", name="Z", lbound=.1, ubound=10),

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
                       
#
# MZ likelihood is set up as pright*(Likelihood|Zygosity=MZ) + pwrong*(Likelihood|DZ)
# DZ likelihood is set up as pright*(Likelihood|Zygosity=DZ) + pwrong*(Likelihood|MZ)
#
# vector=TRUE argument to mxFitFunctionML(),mxExpectationNormal allows mixture distribution on individual likelihoods
#
	mxModel("MZlike",
		mxData(DataMZ, type="raw"), 
		mxExpectationNormal("twinACE.expCovMZ", "twinACE.expMean",c("T1","T2")),
		mxFitFunctionML(vector=T)
	),
	mxModel("DZlike",
		mxData(DataMZ, type="raw"), 
		mxExpectationNormal("twinACE.expCovDZ", "twinACE.expMean",c("T1","T2")),
		mxFitFunctionML(vector=T)
	),
	mxMatrix(type="Full", nrow=dim(frameMZ)[1], ncol=1, values = frameMZ$pMZ, name="pMZ"),
	mxMatrix(type="Unit", nrow=dim(frameMZ)[1], ncol=1, name="unit"),

	mxAlgebra(-2*sum(log(pMZ * MZlike.objective + (unit-pMZ) * DZlike.objective)), name="twin"), 
	mxFitFunctionAlgebra("twin")
)

twinACEFit <- mxRun(twinACEModel)

omxCheckCloseEnough(twinACEFit$output$fit, 4715.549, 1e-2)
