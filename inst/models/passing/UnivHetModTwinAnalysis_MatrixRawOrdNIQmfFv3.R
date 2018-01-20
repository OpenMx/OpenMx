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


# -----------------------------------------------------------------------[]
# Program: UnivHetModTwinAnalysis_MatrixRawCon.R  
#  Author: Hermine Maes
#    Date: 12 01 2009 
#
# Univariate Heterogeneity Twin Analysis model to estimate causes of variation (ACE)
# Matrix style model input - Raw data input
#
# Revision History
#   Hermine Maes -- 02 24 2010 updated & reformatted
#   Hermine Maes -- 02 20 2010 added qualitative differences
# -----------------------------------------------------------------------

require(OpenMx)

# Prepare Data and Print Summary Statistics
# -----------------------------------------------------------------------

stageVars <- c('famno','zyg','zyg2','sex1','sex2','age1','age2','agesex1','agesex2',
'init1','inis1','inits1','init2','inis2','inits2')
sttob <- suppressWarnings(try(read.table("models/passing/data/nicsGarbage.ord",header=F, na.strings=".",col.names=stageVars),
			      silent=TRUE))
if (is(sttob,"try-error")) sttob <- read.table("data/nicsGarbage.ord",header=F, na.strings=".",col.names=stageVars)

Vars <- c('init')
nv <- 1
selVars <- c('init1','init2')   #paste("t",c(rep(1,nv),rep(2,nv)),Vars,sep="")
ntv <- nv*2
mzfDataBin <- subset(sttob, zyg==1, c(selVars,'age1','age2','famno','zyg'))
dzfDataBin <- subset(sttob, zyg==2, c(selVars,'age1','age2','famno','zyg'))
mzmDataBin <- subset(sttob, zyg==3, c(selVars,'age1','age2','famno','zyg'))
dzmDataBin <- subset(sttob, zyg==4, c(selVars,'age1','age2','famno','zyg'))
dzmfDataBin <- subset(sttob, zyg==5, c(selVars,'age1','age2','famno','zyg'))

mzfDataBin[,1:2] <- mxFactor(mzfDataBin[,1:2], levels = c(0:1))  #c(1,2)
dzfDataBin[,1:2] <- mxFactor(dzfDataBin[,1:2], levels = c(0:1))
mzmDataBin[,1:2] <- mxFactor(mzmDataBin[,1:2], levels = c(0:1))
dzmDataBin[,1:2] <- mxFactor(dzmDataBin[,1:2], levels = c(0:1))
dzmfDataBin[,1:2] <- mxFactor(dzmfDataBin[,1:2], levels = c(0:1))
#removeMz <- which(apply(is.na(mzDataBin), 1, sum) == 4)
#removeDz <- which(apply(is.na(dzDataBin), 1, sum) == 4)
#mzDataBin <- mzDataBin[-removeMz,]
#zDataBin <- dzDataBin[-removeDz,]



# Fit Heterogeneity ACE Model with RawData and Matrices Input
# -----------------------------------------------------------------------
univHetACEModelNIQmfFv <- mxModel("univHetACE",
    mxModel("ACE",
    # Matrices a, c, and e to store a, c, and e path coefficients
        mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, values=0.69564726, lbound=.001, ubound=1, label=c("am11"), name="am" ),
        mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, values=0.00445831, lbound=.001, ubound=1, label=c("cm11"), name="cm" ),
        mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, values=0.71836969, lbound=.1, ubound=1, label=c("em11"), name="em" ),
        mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, values=0.73973780, lbound=.001, ubound=1, label=c("af11"), name="af" ),
        mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, values=0.40009607, lbound=.001, ubound=1, label=c("cf11"), name="cf" ),
        mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, values=0.54102784, lbound=.1, ubound=1, label=c("ef11"), name="ef" ),
        mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, values=1, lbound=0, ubound=1, label="rg11", name="rg"),
    # Matrices aI, cI, and eI to store moderated a, c, and e path coefficients
        mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, values=0, label=c("aIm11"), name="aIm" ),
        mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, values=0, label=c("cIm11"), name="cIm" ),
        mxMatrix( type="Lower", nrow=nv, ncol=nv, free=F, values=0, label=c("eIm11"), name="eIm" ),
        mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, values=0, label=c("aIf11"), name="aIf" ),
        mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, values=0, label=c("cIf11"), name="cIf" ),
        mxMatrix( type="Lower", nrow=nv, ncol=nv, free=F, values=0, label=c("eIf11"), name="eIf" ),
    # Matrices A, C, and E compute variance components
        mxAlgebra( expression=am %*% t(am), name="Am" ),
        mxAlgebra( expression=cm %*% t(cm), name="Cm" ),
        mxAlgebra( expression=em %*% t(em), name="Em" ),
        mxAlgebra( expression=af %*% t(af), name="Af" ),
        mxAlgebra( expression=cf %*% t(cf), name="Cf" ),
        mxAlgebra( expression=ef %*% t(ef), name="Ef" ),
    # Algebra to compute total variances and standard deviations (diagonal only)
        mxAlgebra( expression=Am+Cm+Em, name="Vm" ),
        mxAlgebra( expression=Af+Cf+Ef, name="Vf" ),
        mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I"),
        mxAlgebra( expression=solve(sqrt(I*Vm)), name="iSDm"),
        mxAlgebra( expression=solve(sqrt(I*Vf)), name="iSDf"),
    # Constraint on variance of Bininal variables
        mxConstraint( expression=I*Vm==I, name="Var1m"),
        mxConstraint( expression=I*Vf==I, name="Var1f"),
    # Matrix & Algebra for expected means vector and expected thresholds
        mxMatrix( type="Zero", nrow=1, ncol=nv, name="Mean" ),
        mxAlgebra( expression= cbind(Mean,Mean), name="expMean" ),
        mxMatrix( type="Full", nrow=1, ncol=nv, free=FALSE, values=-0.34513213, label=c("threm"), name="Threm" ),
        mxMatrix( type="Full", nrow=1, ncol=nv, free=FALSE, values=-0.26268730, label=c("thref"), name="Thref" ),
        mxMatrix( type="Full", nrow=nv, ncol=2, free=c(FALSE,F), values=c(0.29685484,0), label=c("lm11","qm11"), name="bm" ),
        mxMatrix( type="Full", nrow=nv, ncol=2, free=c(FALSE,F), values=c(-0.19190645,0), label=c("lf11","qf11"), name="bf" ),
    # Matrices & Algebra to plot means by age
        mxMatrix( type="Full", nrow=5, ncol=1, values=c(.3,.4,.5,.6,.7), name="Ages"),
        mxMatrix( type="Full", nrow=2, ncol=5, values=c(.3,.09,.4,.16,.5,.25,.6,.36,.7,.49), name="Agelq"),
        mxMatrix( type="Unit", nrow=5, ncol=1, name="unit"),
        mxAlgebra( expression=unit%x%Threm+ t(bm%*%Agelq), name="TIm"),
        mxAlgebra( expression=unit%x%Thref+ t(bf%*%Agelq), name="TIf"),
        mxAlgebra( expression=rbind(t(vech((am+ 0.3%x%aIm) %*% t(am+ 0.3%x%aIm))),
                                    t(vech((am+ 0.7%x%aIm) %*% t(am+ 0.7%x%aIm)))), name="AIm" ),
        mxAlgebra( expression=rbind(t(vech((cm+ 0.3%x%cIm) %*% t(cm+ 0.3%x%cIm))),
                                    t(vech((cm+ 0.7%x%cIm) %*% t(cm+ 0.7%x%cIm)))), name="CIm" ),
        mxAlgebra( expression=rbind(t(vech((em+ 0.3%x%eIm) %*% t(em+ 0.3%x%eIm))),
                                    t(vech((em+ 0.7%x%eIm) %*% t(em+ 0.7%x%eIm)))), name="EIm" ),
        mxAlgebra( expression=rbind(t(vech((af+ 0.3%x%aIf) %*% t(af+ 0.3%x%aIf))),
                                    t(vech((af+ 0.7%x%aIf) %*% t(af+ 0.7%x%aIf)))), name="AIf" ),
        mxAlgebra( expression=rbind(t(vech((cf+ 0.3%x%cIf) %*% t(cf+ 0.3%x%cIf))),
                                    t(vech((cf+ 0.7%x%cIf) %*% t(cf+ 0.7%x%cIf)))), name="CIf" ),
        mxAlgebra( expression=rbind(t(vech((ef+ 0.3%x%eIf) %*% t(ef+ 0.3%x%eIf))),
                                    t(vech((ef+ 0.7%x%eIf) %*% t(ef+ 0.7%x%eIf)))), name="EIf" ),
        mxAlgebra( expression=AIm+CIm+EIm, name="VIm" ),
        mxAlgebra( expression=cbind(AIm,CIm,EIm,VIm), name="estByAge_Male" ),
        mxAlgebra( expression=cbind(AIm/VIm,CIm/VIm,EIm/VIm,VIm), name="stEstByAge_Male" ),
        mxAlgebra( expression=AIf+CIf+EIf, name="VIf" ),
        mxAlgebra( expression=cbind(AIf,CIf,EIf,VIf), name="estByAge_Female" ),
        mxAlgebra( expression=cbind(AIf/VIf,CIf/VIf,EIf/VIf,VIf), name="stEstByAge_Female" )
    ),
    mxModel("MZm",
    # Matrix for moderating/interacting variable
        mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.age1"), name="Age1"), 
        mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.age2"), name="Age2"), 
    # Matrices A, C, and E compute variance components
        mxAlgebra( expression=(ACE.am+ Age1%x%ACE.aIm) %*% t(ACE.am+ Age1%x%ACE.aIm), name="Am1" ),
        mxAlgebra( expression=(ACE.cm+ Age1%x%ACE.cIm) %*% t(ACE.cm+ Age1%x%ACE.cIm), name="Cm1" ),
        mxAlgebra( expression=(ACE.em+ Age1%x%ACE.eIm) %*% t(ACE.em+ Age1%x%ACE.eIm), name="Em1" ),
        mxAlgebra( expression=(ACE.am+ Age1%x%ACE.aIm) %*% t(ACE.am+ Age2%x%ACE.aIm), name="Am12" ),
        mxAlgebra( expression=(ACE.cm+ Age1%x%ACE.cIm) %*% t(ACE.cm+ Age2%x%ACE.cIm), name="Cm12" ),
        mxAlgebra( expression=(ACE.am+ Age2%x%ACE.aIm) %*% t(ACE.am+ Age1%x%ACE.aIm), name="Am21" ),
        mxAlgebra( expression=(ACE.cm+ Age2%x%ACE.cIm) %*% t(ACE.cm+ Age1%x%ACE.cIm), name="Cm21" ),
        mxAlgebra( expression=(ACE.am+ Age2%x%ACE.aIm) %*% t(ACE.am+ Age2%x%ACE.aIm), name="Am2" ),
        mxAlgebra( expression=(ACE.cm+ Age2%x%ACE.cIm) %*% t(ACE.cm+ Age2%x%ACE.cIm), name="Cm2" ),
        mxAlgebra( expression=(ACE.em+ Age2%x%ACE.eIm) %*% t(ACE.em+ Age2%x%ACE.eIm), name="Em2" ),
    # Algebra for expected variance/covariance matrix and expected mean vector in MZ
        mxAlgebra( expression= rbind  ( cbind(Am1+Cm1+Em1 , Am12+Cm12),
                                        cbind(Am21+Cm21   , Am2+Cm2+Em2)), name="expCovMZm" ),
        mxAlgebra( expression= ACE.bm%*%rbind(Age1,Age1^2), name="AgeRm1"),
        mxAlgebra( expression= ACE.bm%*%rbind(Age2,Age2^2), name="AgeRm2"),
        mxAlgebra( expression= cbind((ACE.Threm + t(AgeRm1)),(ACE.Threm + t(AgeRm2))), dimnames=list('th1',selVars), name="expThreMZm"),
    # Data & Objective
        mxData( observed=mzmDataBin, type="raw" ),
	    mxExpectationNormal( covariance="expCovMZm", means="ACE.expMean", dimnames=selVars, thresholds="expThreMZm"),
	    mxFitFunctionML(vector=TRUE)
    ),
    mxModel("DZm", 
    # Matrix for moderating/interacting variable
        mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.age1"), name="Age1"), 
        mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.age2"), name="Age2"), 
    # Matrices A, C, and E compute variance components        
        mxAlgebra( expression=(ACE.am+ Age1%x%ACE.aIm) %*% t(ACE.am+ Age1%x%ACE.aIm), name="Am1" ),
        mxAlgebra( expression=(ACE.cm+ Age1%x%ACE.cIm) %*% t(ACE.cm+ Age1%x%ACE.cIm), name="Cm1" ),
        mxAlgebra( expression=(ACE.em+ Age1%x%ACE.eIm) %*% t(ACE.em+ Age1%x%ACE.eIm), name="Em1" ),
        mxAlgebra( expression=(ACE.am+ Age1%x%ACE.aIm) %*% t(ACE.am+ Age2%x%ACE.aIm), name="Am12" ),
        mxAlgebra( expression=(ACE.cm+ Age1%x%ACE.cIm) %*% t(ACE.cm+ Age2%x%ACE.cIm), name="Cm12" ),
        mxAlgebra( expression=(ACE.am+ Age2%x%ACE.aIm) %*% t(ACE.am+ Age1%x%ACE.aIm), name="Am21" ),
        mxAlgebra( expression=(ACE.cm+ Age2%x%ACE.cIm) %*% t(ACE.cm+ Age1%x%ACE.cIm), name="Cm21" ),
        mxAlgebra( expression=(ACE.am+ Age2%x%ACE.aIm) %*% t(ACE.am+ Age2%x%ACE.aIm), name="Am2" ),
        mxAlgebra( expression=(ACE.cm+ Age2%x%ACE.cIm) %*% t(ACE.cm+ Age2%x%ACE.cIm), name="Cm2" ),
        mxAlgebra( expression=(ACE.em+ Age2%x%ACE.eIm) %*% t(ACE.em+ Age2%x%ACE.eIm), name="Em2" ),
    # Algebra for expected variance/covariance matrix and expected mean vector in MZ
        mxAlgebra( expression= rbind  ( cbind(Am1+Cm1+Em1     , 0.5%x%Am12+Cm12),
                                        cbind(0.5%x%Am21+Cm21 , Am2+Cm2+Em2)),  name="expCovDZm" ),
        mxAlgebra( expression= ACE.bm%*%rbind(Age1,Age1^2), name="AgeRm1"),
        mxAlgebra( expression= ACE.bm%*%rbind(Age2,Age2^2), name="AgeRm2"),
        mxAlgebra( expression= cbind((ACE.Threm + t(AgeRm1)),(ACE.Threm + t(AgeRm2))), dimnames=list('th1',selVars), name="expThreDZm"),
    # Data & Objective
        mxData( observed=dzmDataBin, type="raw" ),
	    mxExpectationNormal( covariance="expCovDZm", means="ACE.expMean", dimnames=selVars, thresholds="expThreDZm"),
	    mxFitFunctionML(vector=TRUE)
    ),
    mxModel("MZf",
    # Matrix for moderating/interacting variable
        mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.age1"), name="Age1"), 
        mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.age2"), name="Age2"), 
    # Matrices A, C, and E compute variance components
        mxAlgebra( expression=(ACE.af+ Age1%x%ACE.aIf) %*% t(ACE.af+ Age1%x%ACE.aIf), name="Af1" ),
        mxAlgebra( expression=(ACE.cf+ Age1%x%ACE.cIf) %*% t(ACE.cf+ Age1%x%ACE.cIf), name="Cf1" ),
        mxAlgebra( expression=(ACE.ef+ Age1%x%ACE.eIf) %*% t(ACE.ef+ Age1%x%ACE.eIf), name="Ef1" ),
        mxAlgebra( expression=(ACE.af+ Age1%x%ACE.aIf) %*% t(ACE.af+ Age2%x%ACE.aIf), name="Af12" ),
        mxAlgebra( expression=(ACE.cf+ Age1%x%ACE.cIf) %*% t(ACE.cf+ Age2%x%ACE.cIf), name="Cf12" ),
        mxAlgebra( expression=(ACE.af+ Age2%x%ACE.aIf) %*% t(ACE.af+ Age1%x%ACE.aIf), name="Af21" ),
        mxAlgebra( expression=(ACE.cf+ Age2%x%ACE.cIf) %*% t(ACE.cf+ Age1%x%ACE.cIf), name="Cf21" ),
        mxAlgebra( expression=(ACE.af+ Age2%x%ACE.aIf) %*% t(ACE.af+ Age2%x%ACE.aIf), name="Af2" ),
        mxAlgebra( expression=(ACE.cf+ Age2%x%ACE.cIf) %*% t(ACE.cf+ Age2%x%ACE.cIf), name="Cf2" ),
        mxAlgebra( expression=(ACE.ef+ Age2%x%ACE.eIf) %*% t(ACE.ef+ Age2%x%ACE.eIf), name="Ef2" ),
    # Algebra for expected variance/covariance matrix and expected mean vector in MZ
        mxAlgebra( expression= rbind  ( cbind(Af1+Cf1+Ef1 , Af12+Cf12),
                                        cbind(Af21+Cf21   , Af2+Cf2+Ef2)), name="expCovMZf" ),
        mxAlgebra( expression= ACE.bf%*%rbind(Age1,Age1^2), name="AgeRf1"),
        mxAlgebra( expression= ACE.bf%*%rbind(Age2,Age2^2), name="AgeRf2"),
        mxAlgebra( expression= cbind((ACE.Thref + t(AgeRf1)),(ACE.Thref + t(AgeRf2))), dimnames=list('th1',selVars), name="expThreMZf"),
    # Data & Objective
        mxData( observed=mzfDataBin, type="raw" ),
	    mxExpectationNormal( covariance="expCovMZf", means="ACE.expMean", dimnames=selVars, thresholds="expThreMZf"),
	    mxFitFunctionML(vector=TRUE)
    ),
    mxModel("DZf", 
    # Matrix for moderating/interacting variable
        mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.age1"), name="Age1"), 
        mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.age2"), name="Age2"), 
    # Matrices A, C, and E compute variance components        
        mxAlgebra( expression=(ACE.af+ Age1%x%ACE.aIf) %*% t(ACE.af+ Age1%x%ACE.aIf), name="Af1" ),
        mxAlgebra( expression=(ACE.cf+ Age1%x%ACE.cIf) %*% t(ACE.cf+ Age1%x%ACE.cIf), name="Cf1" ),
        mxAlgebra( expression=(ACE.ef+ Age1%x%ACE.eIf) %*% t(ACE.ef+ Age1%x%ACE.eIf), name="Ef1" ),
        mxAlgebra( expression=(ACE.af+ Age1%x%ACE.aIf) %*% t(ACE.af+ Age2%x%ACE.aIf), name="Af12" ),
        mxAlgebra( expression=(ACE.cf+ Age1%x%ACE.cIf) %*% t(ACE.cf+ Age2%x%ACE.cIf), name="Cf12" ),
        mxAlgebra( expression=(ACE.af+ Age2%x%ACE.aIf) %*% t(ACE.af+ Age1%x%ACE.aIf), name="Af21" ),
        mxAlgebra( expression=(ACE.cf+ Age2%x%ACE.cIf) %*% t(ACE.cf+ Age1%x%ACE.cIf), name="Cf21" ),
        mxAlgebra( expression=(ACE.af+ Age2%x%ACE.aIf) %*% t(ACE.af+ Age2%x%ACE.aIf), name="Af2" ),
        mxAlgebra( expression=(ACE.cf+ Age2%x%ACE.cIf) %*% t(ACE.cf+ Age2%x%ACE.cIf), name="Cf2" ),
        mxAlgebra( expression=(ACE.ef+ Age2%x%ACE.eIf) %*% t(ACE.ef+ Age2%x%ACE.eIf), name="Ef2" ),
    # Algebra for expected variance/covariance matrix and expected mean vector in MZ
        mxAlgebra( expression= rbind  ( cbind(Af1+Cf1+Ef1     , 0.5%x%Af12+Cf12),
                                        cbind(0.5%x%Af21+Cf21 , Af2+Cf2+Ef2)),  name="expCovDZf" ),
        mxAlgebra( expression= ACE.bf%*%rbind(Age1,Age1^2), name="AgeRf1"),
        mxAlgebra( expression= ACE.bf%*%rbind(Age2,Age2^2), name="AgeRf2"),
        mxAlgebra( expression= cbind((ACE.Thref + t(AgeRf1)),(ACE.Thref + t(AgeRf2))), dimnames=list('th1',selVars), name="expThreDZf"),
    # Data & Objective
        mxData( observed=dzfDataBin, type="raw" ),
	    mxExpectationNormal( covariance="expCovDZf", means="ACE.expMean", dimnames=selVars, thresholds="expThreDZf"),
	    mxFitFunctionML(vector=TRUE)
    ),

	    mxModel("DZmf", 
    # Matrix for moderating/interacting variable
        mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.age1"), name="Age1"), 
        mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.age2"), name="Age2"), 
    # Matrices A, C, and E compute variance components
        mxAlgebra( expression=(ACE.am+ Age1%x%ACE.aIm) %*% t(ACE.am+ Age1%x%ACE.aIm), name="Am1" ),
        mxAlgebra( expression=(ACE.cm+ Age1%x%ACE.cIm) %*% t(ACE.cm+ Age1%x%ACE.cIm), name="Cm1" ),
        mxAlgebra( expression=(ACE.em+ Age1%x%ACE.eIm) %*% t(ACE.em+ Age1%x%ACE.eIm), name="Em1" ),
        mxAlgebra( expression=(ACE.af+ Age2%x%ACE.aIf) %*% t(ACE.af+ Age2%x%ACE.aIf), name="Af2" ),
        mxAlgebra( expression=(ACE.cf+ Age2%x%ACE.cIf) %*% t(ACE.cf+ Age2%x%ACE.cIf), name="Cf2" ),
        mxAlgebra( expression=(ACE.ef+ Age2%x%ACE.eIf) %*% t(ACE.ef+ Age2%x%ACE.eIf), name="Ef2" ),
        mxAlgebra( expression=(ACE.rg)%x%(ACE.af+ Age2%x%ACE.aIf) %*% t(ACE.am+ Age1%x%ACE.aIm), name="afam" ),
        mxAlgebra( expression=(ACE.cf+ Age2%x%ACE.cIf) %*% t(ACE.cm+ Age1%x%ACE.cIm), name="cfcm" ),
        mxAlgebra( expression=(ACE.rg)%x%(ACE.am+ Age1%x%ACE.aIm) %*% t(ACE.af+ Age2%x%ACE.aIf), name="amaf" ),
        mxAlgebra( expression=(ACE.cm+ Age1%x%ACE.cIm) %*% t(ACE.cf+ Age2%x%ACE.cIf), name="cmcf" ),
    # Algebra for expected variance/covariance matrix and expected mean vector in MZ
        mxAlgebra( expression= rbind  ( cbind(Am1+Cm1+Em1        , 0.5%x%amaf+cmcf),
                                        cbind(0.5%x%afam+cfcm , Af2+Cf2+Ef2)),  name="expCovDZmf" ),
        mxAlgebra( expression= ACE.bm%*%rbind(Age1,Age1^2), name="AgeRm1"),
        mxAlgebra( expression= ACE.bf%*%rbind(Age2,Age2^2), name="AgeRf2"),
        mxAlgebra( expression= cbind((ACE.Threm + t(AgeRm1)),(ACE.Thref + t(AgeRf2))), dimnames=list('th1',selVars), name="expThreDZmf"),
    # Data & Objective
        mxData( observed=dzmfDataBin, type="raw" ),
        mxFitFunctionML(),mxExpectationNormal( covariance="expCovDZmf", means="ACE.expMean", dimnames=selVars, thresholds="expThreDZmf" )
    ),
#univModHetACEFit <- mxModel(ACE,MZm,DZm,MZf,DZf,DZmf,
    mxAlgebra( expression=-2*sum(log(rbind(MZm.objective,DZm.objective,MZf.objective,DZf.objective))), 
        name="minus2sumloglikelihood" ),
    mxFitFunctionAlgebra("minus2sumloglikelihood")
)

univHetACEModelNIQmfFv <- mxOption(univHetACEModelNIQmfFv, "Function precision", 1e-10) 
univHetACEFitNIQmfFv <- mxRun(univHetACEModelNIQmfFv)
omxCheckCloseEnough(univHetACEFitNIQmfFv$output$minimum, 
	univHetACEFitNIQmfFv$output$Minus2LogLikelihood, 0.0001)
