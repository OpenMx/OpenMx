#
#   Copyright 2007-2010 The OpenMx Project
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

DataMZ<-read.table("data/sim1.mz",header<-F)
DataDZ<-read.table("data/sim1.dz",header<-F)
selVars<-c("T1","T2")
names(DataMZ)<-selVars
names(DataDZ)<-selVars
pright<-as.vector(.95)
pwrong<-as.vector(.05)

twinACEModel <- mxModel("twinACE", 
                # Matrix expMean for expected mean vector for MZ and DZ twins    
        mxMatrix("Full", 1, 2, T, 20, "mean", name="expMean"), 

                # Matrices X, Y, and Z to store the a, c, and e path coefficients
        mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=.6, lbound=0, label="a", name="X"),
        mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=.6, lbound=0, label="c", name="Y"),
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
                       
#
# MZ likelihood is set up as pright*(Likelihood|Zygosity=MZ) + pwrong*(Likelihood|DZ)
# DZ likelihood is set up as pright*(Likelihood|Zygosity=DZ) + pwrong*(Likelihood|MZ)
#
# vector=TRUE argument to mxFIMLObjective allows mixture distribution on individual likelihoods
#
        mxModel("MZcorrect",
                mxData(DataMZ, type="raw"), 
                mxFIMLObjective("twinACE.expCovMZ", "twinACE.expMean",selVars, vector=T)),
        mxModel("MZincorrect",
                mxData(DataMZ, type="raw"), 
                mxFIMLObjective("twinACE.expCovDZ", "twinACE.expMean",selVars, vector=T)),

        mxModel("DZcorrect", 
                mxData(DataDZ, type="raw"), 
                mxFIMLObjective("twinACE.expCovDZ", "twinACE.expMean",selVars, vector=T)),
        mxModel("DZincorrect", 
                mxData(DataDZ, type="raw"), 
                mxFIMLObjective("twinACE.expCovMZ", "twinACE.expMean",selVars, vector=T)),


        mxAlgebra(-2*sum(log(pright%x%MZcorrect.objective + pwrong%x%MZincorrect.objective)) + 
                  -2*sum(log(pright%x%DZcorrect.objective + pwrong%x%DZincorrect.objective)), name="twin"), 
        mxAlgebraObjective("twin"))
twinACEFit <- mxRun(twinACEModel)

# Check results against hard-coded Mx1 estimates and likelihood
estimates <- mxEval(as.vector(c(X, Y, Z, expMean[1,1])), twinACEFit)
fitStatistics <- mxEval(objective, twinACEFit)
omxCheckCloseEnough(as.vector(c(0.7974,0.4196,0.4509,-0.0254,10165.966)),as.vector(c(estimates,fitStatistics)),.005)

