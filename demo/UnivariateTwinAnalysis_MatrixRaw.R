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
share <- mxModel("share",
	# Matrices X,Y, and Z to store the a,c,and e path coefficients
    mxMatrix("Full", nrow=1, ncol=1, free=TRUE,  values=.6,  label="a", name="X"), 
    mxMatrix("Full", nrow=1, ncol=1, free=TRUE,  values=.6,  label="c", name="Y"),
    mxMatrix("Full", nrow=1, ncol=1, free=TRUE,  values=.6,  label="e", name="Z"),
    mxAlgebra(X %*% t(X), name="A"), # compute A,C, and E variance components
	mxAlgebra(Y %*% t(Y), name="C"),
	mxAlgebra(Z %*% t(Z), name="E"))

# because both expMeanMZ and expMeanDZ
# share the same label, we are equating them

mzModel <- mxModel(share, name = "MZ",
    mxMatrix("Full", nrow=1, ncol=2, free=TRUE,  values= 20, 
    	label="mean", dimnames=list(NULL, selVars), name="expMean"), 
    # Algebra for expected variance/covariance matrix in MZ
    mxAlgebra(rbind (cbind(A+C+E  , A+C),
                     cbind(A+C    , A+C+E)), 
              dimnames = list(selVars, selVars), name="expCov"),
    mxData(mzfData, type="raw"), 
    mxFIMLObjective("expCov", "expMean"))

dzModel <- mxModel(share, name = "DZ",
    mxMatrix("Full", nrow=1, ncol=2, free=TRUE,  values= 20, 
    	label="mean", dimnames=list(NULL, selVars), name="expMean"),
	# just a constant 0.5 for use in algebras below
    mxMatrix("Full", nrow=1, ncol=1, free=FALSE, values=.5,  name="h"), 
    mxAlgebra(rbind (cbind(A+C+E  , h%x%A+C),
                     cbind(h%x%A+C, A+C+E)), 
              dimnames = list(selVars, selVars), name="expCov"),
    mxData(dzfData, type="raw"), 
    mxFIMLObjective("expCov", "expMean"))
     	

twinACEModel <- mxModel("twinACE", 
	mzModel, dzModel, 
    mxAlgebra(MZ.objective + DZ.objective, name="twin"), 
    mxAlgebraObjective("twin"))

#Run ACE model
twinACEFit <- mxRun(twinACEModel)

MZc <- mxEval(MZ.expCov,  twinACEFit)
DZc <- mxEval(DZ.expCov,  twinACEFit)
M   <- mxEval(MZ.expMean, twinACEFit)

# Retrieve the A, C, and E variance components
A   <- mxEval(MZ.A, twinACEFit)
C   <- mxEval(MZ.C, twinACEFit)
E   <- mxEval(MZ.E, twinACEFit)

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
mzModel <- mxModel(mzModel,
	mxMatrix("Full", nrow=1, ncol=1, free=F, values=0, label="c", name="Y"))

dzModel <- mxModel(dzModel,
	mxMatrix("Full", nrow=1, ncol=1, free=F, values=0, label="c", name="Y"))

twinAEModel <- mxModel("twinAE", 
	mzModel, dzModel, 
    mxAlgebra(MZ.objective + DZ.objective, name="twin"), 
    mxAlgebraObjective("twin"))
	
twinAEFit <- mxRun(twinAEModel)

# As above, retrieve covariance, means, and A, C, and E variance components
MZc <- mxEval(MZ.expCov, twinAEFit)
DZc <- mxEval(DZ.expCov, twinAEFit)
A   <- mxEval(MZ.A, twinAEFit)
C   <- mxEval(MZ.C, twinAEFit)
E   <- mxEval(MZ.E, twinAEFit)
V <- (A + C + E)
a2  <- A / V
c2  <- C / V
e2  <- E / V
# As an example of how R can be used to report information based on OpenMx output, we'll build a reporting table with labels
AEest <- round(rbind(cbind(A, C, E),cbind(a2, c2, e2)),3) # assemble the variance and standardized variance compoents into a table
AEest <- data.frame(AEest, row.names=c("Variance Components","Standardized ∂2")) # build a data.frame with row.names
names(AEest)<-c("A", "C", "E") # add column names
LL_AE <- mxEval(objective, twinAEFit)

LRT_ACE_AE <- LL_AE - LL_ACE

#Print relevant output
ACEest
AEest
LRT_ACE_AE
