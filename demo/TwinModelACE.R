dataMZ <- mxData(matrix(c(1,.8,.8,1), nrow = 2, ncol=2), type="cov", numObs=100)
dataDZ <- mxData(matrix(c(1,.5,.5,1), nrow = 2, ncol=2), type="cov", numObs=100)

matrixX <- mxMatrix("Full",.6,free=TRUE, nrow=1, ncol=1, name="X")
matrixY <- mxMatrix("Full",.6,free=TRUE, nrow=1, ncol=1, name="Y")
matrixZ <- mxMatrix("Full",.6,free=TRUE, nrow=1, ncol=1, name="Z")
matrixh <- mxMatrix("Full",.5,free=FALSE, nrow=1, ncol=1, name="h")
algebraA <- mxAlgebra(X * t(X), name="A")
algebraC <- mxAlgebra(Y * t(Y), name="C")
algebraE <- mxAlgebra(Z * t(Z), name="E")
algebraAC <- mxAlgebra(A + C, name = "AC")
algebraACE <- mxAlgebra(A + C + E, name = "ACE")
algebrahAC <- mxAlgebra(h %x% A + C, name = "hAC")

# Add the shared matrices and algebras to the shared model
sharedModel <- mxModel("share", matrixX, matrixY, matrixZ,
        matrixh, algebraA, algebraC, algebraE,
        algebraAC, algebraACE, algebrahAC)

covarianceMZ <- mxAlgebra(rbind(
       cbind(share.ACE, share.AC),
       cbind(share.AC, share.ACE)), name="cMZ")

covarianceDZ <- mxAlgebra(rbind(
       cbind(share.ACE, share.hAC),
       cbind(share.hAC, share.ACE)), name="cDZ")

objMZ <- mxMLObjective(covariance = "cMZ")
objDZ <- mxMLObjective(covariance = "cDZ")

modelMZ <- mxModel("modelMZ", dataMZ, covarianceMZ, objMZ)
modelDZ <- mxModel("modelDZ", dataDZ, covarianceDZ, objDZ)

twin <- mxAlgebra(modelMZ.objective + modelDZ.objective, name="twin")
obj <- mxAlgebraObjective("twin")

# Finally add the submodels to the shared model
sharedModel <- mxModel(sharedModel, twin, obj, modelMZ, modelDZ)

# Run the model
sharedModelOut <- mxRun(sharedModel)

expectedACE <- c(.6, .2, .2)
observedACE <- mxEvaluate(c(A, C, E), sharedModelOut)

omxCheckCloseEnough(expectedACE, observedACE, epsilon = 10 ^ -4)

