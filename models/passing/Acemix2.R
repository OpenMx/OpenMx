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
selVars<-c("T1","T2","pMZ")
names(DataMZ)<-selVars
frameMZ<-data.frame(pMZ<-DataMZ$pMZ)

twinACEModel <- mxModel("twinACE", 
                # Matrix expMean for expected mean vector for MZ and DZ twins    
        mxMatrix("Full", 1, 2, T, 20, "mean", name="expMean"), 

                # Matrices X, Y, and Z to store the a, c, and e path coefficients
        mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=.6, label="a", name="X"),
        mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=.6, label="c", name="Y"),
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
        mxModel("MZlike",
                mxData(DataMZ, type="raw"), 
                mxFIMLObjective("twinACE.expCovMZ", "twinACE.expMean",c("T1","T2"), vector=T)),
        mxModel("DZlike",
                mxData(DataMZ, type="raw"), 
                mxFIMLObjective("twinACE.expCovDZ", "twinACE.expMean",c("T1","T2"), vector=T)),

        mxMatrix(type="Full", nrow=dim(frameMZ)[1], ncol=1, values=frameMZ$pMZ, name="pMZ"),
        mxMatrix(type="Unit", nrow=dim(frameMZ)[1], ncol=1, name="unit"),
        
        mxAlgebra(-2*sum(log(pMZ * MZlike.objective + (unit-pMZ) * DZlike.objective)), name="twin"), 
        mxAlgebraObjective("twin")
        )

twinACEFit <- mxRun(twinACEModel)

