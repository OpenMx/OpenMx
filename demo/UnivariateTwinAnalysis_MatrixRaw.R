require(OpenMx)

#Prepare Data
data(twinData)
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
  mxAlgebra(Z %*% t(Z), name="E"),
  mxMatrix("Full", nrow=1, ncol=2, free=TRUE, values= 20, 
     label="mean", name="expMean")
)

mzModel <- mxModel(name = "MZ",
    # Algebra for expected variance/covariance matrix in MZ
    mxAlgebra(rbind (cbind(share.A + share.C + share.E  ,share.A + share.C),
                     cbind(share.A + share.C    ,share.A + share.C + share.E)), 
              name="expCov"),
    mxData(mzfData, type="raw"), 
    mxFIMLObjective("expCov", "share.expMean", selVars)
)

dzModel <- mxModel(name = "DZ",
    # note use of 0.5, converted to 1*1 matrix
    mxAlgebra(rbind (cbind(share.A + share.C + share.E  , 0.5 %x% share.A + share.C),
                     cbind(0.5 %x% share.A + share.C, share.A + share.C + share.E)), 
              name="expCov"), 
    mxData(dzfData, type="raw"), 
    mxFIMLObjective("expCov", "share.expMean", selVars)
)
     	

twinACEModel <- mxModel("twinACE", 
	share, mzModel, dzModel, 
    mxAlgebra(MZ.objective + DZ.objective, name="twin"), 
    mxAlgebraObjective("twin"))

#Run ACE model
twinACEFit <- mxRun(twinACEModel)

MZc <- mxEval(MZ.expCov,  twinACEFit)
DZc <- mxEval(DZ.expCov,  twinACEFit)
M   <- mxEval(share.expMean, twinACEFit)

# Retrieve the A, C, and E variance components
A   <- mxEval(share.A, twinACEFit)
C   <- mxEval(share.C, twinACEFit)
E   <- mxEval(share.E, twinACEFit)

totalVariance <- (A+C+E)
a2  <- A/totalVariance  # Standardize the variance components
c2  <- C/totalVariance
e2  <- E/totalVariance
# As an example of how R can be used to report information based on OpenMx output, we'll build a reporting table with labels
ACEest <- rbind(cbind(A,C,E),cbind(a2,c2,e2)) # join the variance component matrices into rows, and join the rows into a table
ACEest <- data.frame(ACEest, row.names=c("Variance Components","Standardized ∂2")) # build a data.frame with row.names
names(ACEest)<-c("A", "C", "E") # add column names

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
LL_ACE <- mxEval(objective, twinACEFit) # extract loglikelihood
#Compare OpenMx results to Mx results (LL: likelihood; EC: expected covariance, EM: expected means)
omxCheckCloseEnough(LL_ACE,Mx.LL_ACE,.001)
omxCheckCloseEnough(A,Mx.A,.001)
omxCheckCloseEnough(C,Mx.C,.001)
omxCheckCloseEnough(E,Mx.E,.001)
omxCheckCloseEnough(M,Mx.M,.001)

ACEest; LL_ACE; # print table of results and LL

########## RUN a reduced AE MODEL ###########
share <- mxModel(share, mxMatrix("Full", nrow=1, ncol=1, free=FALSE, values=0, label="c", name="Y") ) # drop c at 0

twinAEModel <- mxModel(twinACEModel, share, name = "twinAE")
	
twinAEFit <- mxRun(twinAEModel)

# As above, retrieve covariance, means, and A, C, and E variance components
MZc <- mxEval(MZ.expCov, twinAEFit)
DZc <- mxEval(DZ.expCov, twinAEFit)
A   <- mxEval(share.A, twinAEFit)
C   <- mxEval(share.C, twinAEFit)
E   <- mxEval(share.E, twinAEFit)
totalVariance <- (A + C + E)
a2  <- A / totalVariance
c2  <- C / totalVariance
e2  <- E / totalVariance
# Build a reporting table with labels
AEest <- round(rbind(cbind(A, C, E),cbind(a2, c2, e2)),3) # assemble the variance and standardized variance compoents into a table
AEest <- data.frame(AEest, row.names=c("Variance Components","Standardized ∂2")) # build a data.frame with row.names
names(AEest) <- c("A", "C", "E") # add column names

LL_AE <- mxEval(objective, twinAEFit); # get the log-likelihood
LRT_ACE_AE <- LL_AE - LL_ACE
#Print relevant output
AEest
LRT_ACE_AE
