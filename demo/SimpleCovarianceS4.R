library(OpenMx)

testFit <- function(objective, startVals=c(), bounds=c(), matList=list(), varList=list(), algList=list(), data=c(), state=c()) {
	return(.Call("callNPSOL", objective, startVals, bounds, matList, varList, algList, data, state));
}

model <- mxModel()
model <- mxModel(model, mxMatrix("Full", c(0,0.2,0,0), name = "A", nrow = 2, ncol = 2))
model <- mxModel(model, mxMatrix("Full", c(0.8,0,0,0.8), name="S", nrow=2, ncol=2, free=TRUE))
model <- mxModel(model, mxMatrix("Full", c(1,0,0,1), name="F", nrow=2, ncol=2))
model <- mxModel(model, mxMatrix("Iden", name="I", nrow=2))
covMatrix <- matrix(c(0.77642931, 0.39590663, 0.39590663, 0.49115615),
	nrow = 2, ncol = 2, byrow = TRUE)

model[["A"]]@specification[2,1] <- NA
model[["S"]]@specification[2,1] <- 0
model[["S"]]@specification[1,2] <- 0
model[["S"]]@specification[1,1] <- "apple"
model[["S"]]@specification[2,2] <- "banana"

pList <- omxGenerateParameterList(model)
mList <- omxGenerateSimpleMatrixList(model)
startVals <- omxGenerateValueList(model)
optype <- "cov"
mList[[4]] <- covMatrix
mxCov <- mxMatrix("Symm", covMatrix, nrow=2, ncol=2)
objective <- mxObjective("RAM", covariance = mxCov)
class(objective)
NPSOLOutput <- testFit(objective, startVals, matList = mList, varList = pList, data=covMatrix)
#NPSOLOutput <- .Call("callNPSOL", optype, startVals, list(), mList, pList, NA, NA, NA, new.env())
print(NPSOLOutput)
