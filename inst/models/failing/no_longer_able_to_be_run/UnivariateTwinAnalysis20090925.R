# -----------------------------------------------------------------------
# Program: UnivariateTwinAnalysis20090925.R  
#  Author: Hermine Maes
#    Date: Wed Sep 25 11:45:52 EDT 2009
#
# Revision History
#   Hermine Maes -- Wed Sep 25 11:45:52 EDT 2009 UnivariateTwinAnalysis20090925.R
# 	TBATES: 2017-04-14 04:18PM
# 	Fixed script bug (DataMZ instead of MZ and DataDZ instead of DZ as the data names)
# 	Replace out of date mxFIMLObjective calls
# 	add mxFitFunctionMultigroup(c("MZ", "DZ"))
#   Gave up as there are additional script errors (like using the model object as a string)
# 	and no erorr is identified in this script.
# TODO include equivalent file for mx 1.x 
# TODO: Add omxCheckCloseEnough() calls
# -----------------------------------------------------------------------


# Simulate Data: two standardized variables t1 & t2 for MZ's & DZ's
# -----------------------------------------------------------------------
require(OpenMx)
require(MASS)

set.seed(200)
a2<-0.5		#Additive genetic variance component (a squared)
c2<-0.3		#Common environment variance component (c squared)
e2<-0.2		#Specific environment variance component (e squared)
rMZ <- a2+c2
rDZ <- .5*a2+c2
DataMZ <- mvrnorm(1000, c(0,0), matrix(c(1,rMZ,rMZ,1),2,2))
DataDZ <- mvrnorm(1000, c(0,0), matrix(c(1,rDZ,rDZ,1),2,2))

selVars <- c('t1','t2')
dimnames(DataMZ) <- list(NULL,selVars)
dimnames(DataDZ) <- list(NULL,selVars)
summary(DataMZ)
summary(DataDZ)
colMeans(DataMZ, na.rm=TRUE)
colMeans(DataDZ, na.rm=TRUE)
cov(DataMZ,use="complete")
cov(DataDZ,use="complete")

# Specify and Run Saturated Model with RawData and Matrix-style Input
# -----------------------------------------------------------------------
twinSatModel <- mxModel("twinSat",
	mxModel("MZ",
		mxMatrix("Full", 1, 2, T, c(0,0), dimnames=list(NULL, selVars), name="expMeanMZ"), 
		mxMatrix("Lower", 2, 2, T, .5, dimnames=list(selVars, selVars), name="CholMZ"), 
		mxAlgebra(CholMZ %*% t(CholMZ), name="expCovMZ", dimnames=list(selVars, selVars)), 
		mxData(DataMZ, type="raw"),
		mxExpectationNormal("expCovMZ", "expMeanMZ"),
		mxFitFunctionML()

	),  
	mxModel("DZ",
		mxMatrix("Full", 1, 2, T, c(0,0), dimnames=list(NULL, selVars), name="expMeanDZ"), 
		mxMatrix("Lower", 2, 2, T, .5, dimnames=list(selVars, selVars), name="CholDZ"), 
		mxAlgebra(CholDZ %*% t(CholDZ), name="expCovDZ", dimnames=list(selVars, selVars)), 
		mxData(DataDZ, type="raw"), 
		mxExpectationNormal("expCovDZ", "expMeanDZ"),
		mxFitFunctionML()
	),
	mxFitFunctionMultigroup(c("MZ", "DZ"))
)
twinSatFit <- mxRun(twinSatModel)

# Generate Saturated Model Output
# -----------------------------------------------------------------------
ExpMeanMZ <- mxEval(MZ.expMeanMZ, twinSatFit)
ExpCovMZ  <- mxEval(MZ.expCovMZ,  twinSatFit)
ExpMeanDZ <- mxEval(DZ.expMeanDZ, twinSatFit)
ExpCovDZ  <- mxEval(DZ.expCovDZ,  twinSatFit)
LL_Sat    <- mxEval(objective,    twinSatFit)


# Specify and Run Saturated SubModel 1 equating means across twin order
# -----------------------------------------------------------------------
twinSatModelSub1 <- mxModel(twinSatModel,
	mxModel(MZ, mxMatrix("Full", 1, 2, T, 0, "mMZ", dimnames=list(NULL, selVars), name = "expMeanMZ")), 
	mxModel(DZ, mxMatrix("Full", 1, 2, T, 0, "mDZ", dimnames=list(NULL, selVars), name = "expMeanDZ"))
)
twinSatFitSub1 <- mxRun(twinSatModelSub1)

# Specify and Run Saturated SubModel 2 equating means across twin order and zygosity
# -----------------------------------------------------------------------
twinSatModelSub2 <- mxModel(twinSatModelSub1,
	mxModel("MZ",
		mxMatrix("Full", 1, 2, T, 0, "mean", dimnames=list(NULL, selVars), name="expMeanMZ"), 
		mxMatrix("Lower", 2, 2, T, .5, labels= c("var","MZcov","var"), 
		    dimnames=list(selVars, selVars), name="CholMZ")
	), 
	mxModel("DZ", 
		mxMatrix("Full", 1, 2, T, 0, "mean", dimnames=list(NULL, selVars), name="expMeanDZ"), 
		mxMatrix("Lower", 2, 2, T, .5, labels= c("var","DZcov","var"), 
		    dimnames=list(selVars, selVars), name="CholDZ")
	)
)
twinSatFitSub2 <- mxRun(twinSatModelSub2)

# Generate Saturated Model Comparison Output
# -----------------------------------------------------------------------
LL_Sat  <- mxEval(objective, twinSatFit)
LL_Sub1 <- mxEval(objective, twinSatFitSub1)
LRT1    <- LL_Sub1 - LL_Sat
LL_Sub2 <- mxEval(objective, twinSatFitSub1)
LRT2    <- LL_Sub2 - LL_Sat


# Specify and Run ACE Model with RawData and Matrix-style Input
# -----------------------------------------------------------------------
twinACEModel <- mxModel("twinACE", 
	mxMatrix("Full", 1, 2, T, 20, "mean", dimnames=list(NULL, selVars), name="expMean"), 
		# Matrix expMean for expected mean vector for MZ and DZ twins    
	mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=.6, label="a", name="X"),
	mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=.6, label="c", name="Y"),
	mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=.6, label="e", name="Z"),
		# Matrices X, Y, and Z to store the a, c, and e path coefficients
	mxMatrix("Full", nrow=1, ncol=1, free=FALSE, values=.5, name="h"),
	mxAlgebra(X * t(X), name="A"),
	mxAlgebra(Y * t(Y), name="C"),
	mxAlgebra(Z * t(Z), name="E"),
		# Matrixes A, C, and E to compute A, C, and E variance components
	mxAlgebra(rbind(cbind(A+C+E   , A+C),
					cbind(A+C     , A+C+E)), dimnames = list(selVars, selVars), name="expCovMZ"),
		# Matrix expCOVMZ for expected covariance matrix for MZ twins
	mxAlgebra(rbind(cbind(A+C+E   , h%x%A+C),
					cbind(h%x%A+C , A+C+E)), dimnames = list(selVars, selVars), name="expCovDZ"),
		# Matrix expCOVMZ for expected covariance matrix for DZ twins
	mxModel("MZ",
		mxData(DataMZ, type="raw"), 
		mxFIMLObjective("twinACE.expCovMZ", "twinACE.expMean")),
	mxModel("DZ", 
		mxData(DataDZ, type="raw"), 
		mxFIMLObjective("twinACE.expCovDZ", "twinACE.expMean")),
	mxAlgebra(MZ.objective + DZ.objective, name="twin"), 
	mxFitFunctionAlgebra("twin")
)
twinACEFit <- mxRun(twinACEModel)

# Generate ACE Model Output
# -----------------------------------------------------------------------
LL_ACE <- mxEval(objective, twinACEFit)
LRT_ACE= LL_ACE - LL_Sat

#Retrieve expected mean vector and expected covariance matrices
	MZc <- mxEval(expCovMZ, twinACEFit)
	DZc <- mxEval(expCovDZ, twinACEFit)
	M   <- mxEval(expMean, twinACEFit)
#Retrieve the A, C, and E variance components
	A <- mxEval(A, twinACEFit)
	C <- mxEval(C, twinACEFit)
	E <- mxEval(E, twinACEFit)
#Calculate standardized variance components
	V  <- (A+C+E)
	a2 <- A/V
	c2 <- C/V
	e2 <- E/V
#Build and print reporting table with row and column names
	ACEest <- rbind(cbind(A,C,E),cbind(a2,c2,e2)) 
	ACEest <- data.frame(ACEest, row.names=c("Variance Components","Standardized VC"))
	names(ACEest)<-c("A", "C", "E")
 	ACEest; LL_ACE; LRT_ACE

# Specify and reduced AE Model (drop c $0)
# -----------------------------------------------------------------------
twinAEModel <- mxModel(twinACEModel, name="twinAE",
    mxMatrix("Full", nrow=1, ncol=1, free=F, values=0, label="c", name="Y")
)
twinAEFit <- mxRun(twinAEModel)

# Generate ACE Model Output
# -----------------------------------------------------------------------
LL_AE <- mxEval(objective, twinAEFit)
#Retrieve expected mean vector and expected covariance matrices
	MZc <- mxEval(expCovMZ, twinAEFit)
	DZc <- mxEval(expCovDZ, twinAEFit)
	M   <- mxEval(expMean, twinAEFit)
#Retrieve the A, C and E variance components
	A <- mxEval(A, twinAEFit)
	C <- mxEval(C, twinAEFit)
	E <- mxEval(E, twinAEFit)
#Calculate standardized variance components
	V <- (A+C+E)
	a2 <- A/V
	c2 <- C/V
	e2 <- E/V
#Build and print reporting table with row and column names
	AEest <- rbind(cbind(A,C,E),cbind(a2,c2,e2)) 
	AEest <- data.frame(ACEest, row.names=c("Variance Components","Standardized VC"))
	names(ACEest)<-c("A", "C", "E")
	AEest; LL_AE; 
#Calculate and print likelihood ratio test
	LRT_ACE_AE <- LL_AE - LL_ACE
	LRT_ACE_AE
