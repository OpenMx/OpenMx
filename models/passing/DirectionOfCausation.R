# Program: Direction of Causation Twin Models

require(OpenMx)
nv<-2
# Fit Direction of Causation Model with Cell Frequencies Input
# ---------------------------------------------------------------------
DirectionOfCausationModel <- mxModel("DirectionOfCausation",
    mxModel("ACE",
    # Matrices a, c, and e to store a, c, and e path coefficients
        mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=.6,  name="a" ),
        mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=.6,  name="c" ),
        mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=sqrt(.28), name="e" ),
        mxMatrix( type="Full", nrow=nv, ncol=nv, free=c(F,T,F,F), name="B"),
    # Matrices A, C, and E compute variance components
        mxAlgebra( expression=a %*% t(a), name="A" ),
        mxAlgebra( expression=c %*% t(c), name="C" ),
        mxAlgebra( expression=e %*% t(e), name="E" ),
    # Algebra to compute total variances and standard deviations (diagonal only)
        mxAlgebra( expression=diag2vec(A+C+E), name="V" ),
        mxMatrix( type="Unit", nrow=nv, ncol=1, name="U"),
        mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I"),
        mxAlgebra( expression=solve(vec2diag(sqrt(V))), name="sd"),
        mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, values=0, labels=c("thrV1","thrV2"), name="threshold"),
        mxMatrix( type="Full", nrow=1, ncol=nv*2, free=F, values=Inf,name="infinities"),
        mxMatrix( type="Full", nrow=1, ncol=nv*2, free=F, values=-Inf,name="neginfinities"),
        mxAlgebra( expression=rbind(neginfinities,cbind(threshold,threshold),infinities), name="thresholds"),
        
    # Constraint on variance of A+C+E latent variables
        mxConstraint(V == U, name="Var1"),
    # Algebra for expected variance/covariance matrix in MZ
        mxAlgebra( expression= (I%x%solve(I-B)) %&% rbind  ( cbind(A+C+E , A+C),
                                        cbind(A+C   , A+C+E)), name="expCovMZ" ),
    # Algebra for expected variance/covariance matrix in DZ
        mxAlgebra( expression= (I%x%solve(I-B)) %&% rbind  ( cbind(A+C+E     , 0.5%x%A+C),
                                        cbind(0.5%x%A+C , A+C+E)),  name="expCovDZ"),
    # Mean
        mxMatrix( type="Zero", nrow=1, ncol=nv, name="M" ),
        mxAlgebra( expression= cbind(M,M), name="expMean" ),
	# Integrals
        mxAlgebra(omxAllInt(expCovMZ, expMean, thresholds), name="AllintMZ"),
        mxAlgebra(omxAllInt(expCovDZ, expMean, thresholds), name="AllintDZ"),
		mxAlgebra(cbind(AllintMZ[1:4,],AllintMZ[5:8,],AllintMZ[9:12,],AllintMZ[13:16,]),name="MAllintMZ"),
		mxAlgebra(cbind(AllintDZ[1:4,],AllintDZ[5:8,],AllintDZ[9:12,],AllintDZ[13:16,]),name="MAllintDZ"),
		mxAlgebra(vech(MAllintMZ+t(MAllintMZ)-vec2diag(diag2vec(MAllintMZ))),name="MZExpectedFrequencies"),
		mxAlgebra(vech(MAllintDZ+t(MAllintDZ)-vec2diag(diag2vec(MAllintDZ))),name="DZExpectedFrequencies"),
		
    mxModel("MZ",
		mxMatrix(type="Full", nrow=10, ncol=1, free=FALSE, values=c(141,35,32,25,15,7,33,18,39,47), 			name="MZObservedFrequencies"),
        mxAlgebra( -2 * sum(MZObservedFrequencies * log(ACE.MZExpectedFrequencies)),name="MZalgobj"),		mxAlgebraObjective("MZalgobj")),
    mxModel("DZ", 
		mxMatrix(type="Full", nrow=10, ncol=1, free=FALSE, values=c(58,18,27,44,7,6,33,15,38,81), 			name="DZObservedFrequencies"),		
        mxAlgebra( -2 * sum(DZObservedFrequencies * log(ACE.DZExpectedFrequencies)),name="DZalgobj"),		mxAlgebraObjective("DZalgobj"))),
    mxAlgebra( MZ.objective + DZ.objective, name="neg2sumll" ),
    mxAlgebraObjective("neg2sumll"))
     
DirectionOfCausationRun<-mxRun(DirectionOfCausationModel)     
summary(DirectionOfCausationRun)

                                        
                                        
