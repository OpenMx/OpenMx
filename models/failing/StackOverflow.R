require(OpenMx)

dataMZ <- mxData(matrix(c(1,.8,.8,1), nrow = 2, ncol=2), type="cov",
numObs=100)
dataDZ <- mxData(matrix(c(1,.5,.5,1), nrow = 2, ncol=2), type="cov",
numObs=100)

X <- mxMatrix("Full",.6,free=TRUE, labels='X', nrow=1, ncol=1, name="X")
Y <- mxMatrix("Full",.6,free=TRUE, labels='Y', nrow=1, ncol=1, name="Y")
Z <- mxMatrix("Full",.6,free=TRUE, labels='Z', nrow=1, ncol=1, name="Z")
h <- mxMatrix("Full",.5,free=FALSE, nrow=1, ncol=1, name="h")
A <- mxAlgebra(X * t(X), name="A")
C <- mxAlgebra(Y * t(Y), name="C")
E <- mxAlgebra(Z * t(Z), name="E")
cMZ <- mxAlgebra(rbind(cbind(A+C+E,A+C),cbind(A+C,A+C+E)), name="cMZ")
cDZ <- mxAlgebra(rbind(cbind(A+C+E,h%x%A+C),cbind(h%x%A+C,A+C+E)),
name="cDZ")

objMZ <- mxMLObjective("cMZ")
objDZ <- mxMLObjective("cDZ")

modelMZ <- mxModel("modelMZ", dataMZ, X,Y,Z,h,A,C,E,cMZ, objMZ)
modelDZ <- mxModel("modelDZ", dataDZ, X,Y,Z,h,A,C,E,cDZ, objDZ)

twin <- mxAlgebra(modelMZ.objective + modelDZ.objective, name="twin")
obj <- mxAlgebraObjective("twin")

model <- mxModel("both", twin, obj, modelMZ, modelDZ)
model <- mxRun(model)