# This demo is a variant of TwinModelACE that
# takes advantage of algebraic substitution

dataMZ <- mxData(matrix(c(1,.8,.8,1), nrow = 2, ncol=2), type="cov", numObs=100)
dataDZ <- mxData(matrix(c(1,.5,.5,1), nrow = 2, ncol=2), type="cov", numObs=100)

X <- mxMatrix("Full",.6,free=TRUE, nrow=1, ncol=1, name="X")
Y <- mxMatrix("Full",.6,free=TRUE, nrow=1, ncol=1, name="Y")
Z <- mxMatrix("Full",.6,free=TRUE, nrow=1, ncol=1, name="Z")
h <- mxMatrix("Full",.5,free=FALSE, nrow=1, ncol=1, name="h")
A <- mxAlgebra(X * t(X), name="A")
C <- mxAlgebra(Y * t(Y), name="C")
E <- mxAlgebra(Z * t(Z), name="E")
AC <- mxAlgebra(share.A + share.C, name = "AC")
ACE <- mxAlgebra(share.A + share.C + share.E, name = "ACE")
hAC <- mxAlgebra(share.h %x% share.A + share.C, name = "hAC")

# Add the shared matrices and algebras to the shared model
sharedModel <- mxModel("share", X, Y, Z, h, A, C, E, AC, ACE, hAC)

cMZ <- mxMatrix("Symm", free = FALSE, nrow = 2, ncol = 2,
	labels = (c('share.ACE', 'share.AC', 'share.ACE')), name = 'cMZ')

cDZ <- mxMatrix("Symm", free = FALSE, nrow = 2, ncol = 2,
	labels = (c('share.ACE', 'share.hAC', 'share.ACE')), name = 'cDZ')

objMZ <- mxMLObjective(covariance = "cMZ")
objDZ <- mxMLObjective(covariance = "cDZ")

modelMZ <- mxModel("modelMZ", dataMZ, cMZ, objMZ)
modelDZ <- mxModel("modelDZ", dataDZ, cDZ, objDZ)

twin <- mxAlgebra(modelMZ.objective + modelDZ.objective, name="twin")
obj <- mxAlgebraObjective("twin")

# Finally add the submodels to the shared model
sharedModel <- mxModel(sharedModel, twin, obj, modelMZ, modelDZ)

# Run the model
sharedModelOut <- mxRun(sharedModel)

expectedACE <- c(.6, .2, .2)
observedACE <- c(sharedModelOut$A@result, sharedModelOut$C@result, sharedModelOut$E@result)

omxCheckCloseEnough(expectedACE, observedACE, epsilon = 10 ^ -4)