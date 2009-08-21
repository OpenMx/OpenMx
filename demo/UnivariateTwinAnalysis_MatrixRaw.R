require(OpenMx)

#Prepare Data
data(twinData)
twinVars <- c('fam','age','zyg','part','wt1','wt2','ht1','ht2','htwt1','htwt2','bmi1','bmi2')
#dimnames(twinData) <- list(NULL, twinVars)
summary(twinData)
selVars <- c('bmi1','bmi2')
mzfData <- as.matrix(subset(twinData, zyg==1, c(bmi1,bmi2)))
dzfData <- as.matrix(subset(twinData, zyg==3, c(bmi1,bmi2)))
colMeans(mzfData,na.rm=TRUE)
colMeans(dzfData,na.rm=TRUE)
cov(mzfData,use="complete")
cov(dzfData,use="complete")

#Fit ACE Model with RawData and Matrices Input
twinACEModel <- mxModel("twinACE", 
    mxMatrix("Full", nrow=1, ncol=2, free=TRUE,  values= 20, label="mean", dimnames=list(NULL, selVars), name="expMeanMZ"), 
    mxMatrix("Full", nrow=1, ncol=2, free=TRUE,  values= 20, label="mean", dimnames=list(NULL, selVars), name="expMeanDZ"), 
    mxMatrix("Full", nrow=1, ncol=1, free=TRUE,  values=.6,  label="a", name="X"),
    mxMatrix("Full", nrow=1, ncol=1, free=TRUE,  values=.6,  label="c", name="Y"),
    mxMatrix("Full", nrow=1, ncol=1, free=TRUE,  values=.6,  label="e", name="Z"),
    mxMatrix("Full", nrow=1, ncol=1, free=FALSE, values=.5,  name="h"),
    mxAlgebra(X %*% t(X), name="A"),
    mxAlgebra(Y %*% t(Y), name="C"),
    mxAlgebra(Z %*% t(Z), name="E"),
    # Algebra for expected variance/covariance matrix in MZs
    mxAlgebra(rbind (cbind(A+C+E  , A+C),
                     cbind(A+C    , A+C+E)), dimnames = list(selVars, selVars), name="expCovMZ"),
    # Algebra for expected variance/covariance matrix in DZs
    mxAlgebra(rbind (cbind(A+C+E  , h%x%A+C),
                     cbind(h%x%A+C, A+C+E)), dimnames = list(selVars, selVars), name="expCovDZ"),
		# Build model for the MZ data.
    mxModel("MZ",
			  # read in the manifest variables (selVars), to provide the observed covariance matrix and means for this group. 
        mxData(mzfData, type="raw"), 
				# link these observations to an objective - our expected MZ covariance and means, 
				# so that the likelihood of the observed data can be calculated from their departure from our expectations.
        mxFIMLObjective("twinACE.expCovMZ", "twinACE.expMeanMZ")
    ),
		# Build model for the DZ data.
    mxModel("DZ", 
        mxData(dzfData, type="raw"), 
        mxFIMLObjective("twinACE.expCovDZ", "twinACE.expMeanDZ")),
    mxAlgebra(MZ.objective + DZ.objective, name="twin"), 
    mxAlgebraObjective("twin"))

#Run ACE model
twinACEFit <- mxRun(twinACEModel)

MZc <- mxEval(expCovMZ, twinACEFit)
DZc <- mxEval(expCovDZ, twinACEFit)
M   <- mxEval(expMeanMZ, twinACEFit)

# Retrieve the A, C, and E variance components
A   <- mxEval(A, twinACEFit)
C   <- mxEval(C, twinACEFit)
E   <- mxEval(E, twinACEFit)

totalVariance <- (A+C+E)
a2  <- A/totalVariance  # Standardize the variance components
c2  <- C/totalVariance
e2  <- E/totalVariance
ACEest <- rbind(cbind(A,C,E),cbind(a2,c2,e2))
LL_ACE <- mxEval(objective, twinACEFit)


#Run Mx scripts
#mymatrices1 <- omxOriginalMx("mx-scripts/univACE.mx", "temp-files")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
#Mx.A <- mymatrices1$A1.1
#Mx.C <- mymatrices1$C1.1
#Mx.E <- mymatrices1$E1.1
#Mx.M <- mymatrices1$M2.1
#Mx.LL_ACE <- mymatrices1$F4.1;

#Mx answers hard-coded
#1: Heterogeneity Model
Mx.A <- 0.6173023
Mx.C <- 5.595822e-14
Mx.E <- 0.1730462
Mx.M <- matrix(c(21.39293, 21.39293),1,2)
Mx.LL_ACE <- 4067.663


#Compare OpenMx results to Mx results (LL: likelihood; EC: expected covariance, EM: expected means)
omxCheckCloseEnough(LL_ACE,Mx.LL_ACE,.001)
omxCheckCloseEnough(A,Mx.A,.001)
omxCheckCloseEnough(C,Mx.C,.001)
omxCheckCloseEnough(E,Mx.E,.001)
omxCheckCloseEnough(M,Mx.M,.001)


#Run AE model
twinAEModel <- mxModel(twinACEModel, 
    mxMatrix("Full", nrow=1, ncol=1, free=F, values=0, label="c", name="Y")
    )
twinAEFit <- mxRun(twinAEModel)

# As above, retrieve covariance, means, and A, C, and E variance components
MZc <- mxEval(expCovMZ, twinAEFit)
DZc <- mxEval(expCovDZ, twinAEFit)
A   <- mxEval(A, twinAEFit)
C   <- mxEval(C, twinAEFit)
E   <- mxEval(E, twinAEFit)
V <- (A + C + E)
a2  <- A / V
c2  <- C / V
e2  <- E / V
AEest <- rbind(cbind(A, C, E),cbind(a2, c2, e2))
LL_AE <- mxEval(objective, twinAEFit)

LRT_ACE_AE <- LL_AE - LL_ACE

#Print relevant output
ACEest
AEest
LRT_ACE_AE
