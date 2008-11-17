library(OpenMx)

model <- new("MxModel")
A <- mxMatrix("Full", c(0,0.2,0,0), nrow = 2, ncol = 2)
S <- mxMatrix("Full", c(0.8,0,0,0.8), nrow=2, ncol=2, free=TRUE)
F <- mxMatrix("Full", c(1,0,0,1), nrow=2, ncol=2)
I <- mxMatrix("Iden", nrow=2)
covMatrix <- matrix(c(0.77642931, 0.39590663, 0.39590663, 0.49115615),
	nrow = 2, ncol = 2, byrow = TRUE)

A@specification[2,1] <- NA
S@specification[2,1] <- 0
S@specification[1,2] <- 0
S@specification[1,1] <- "apple"
S@specification[2,2] <- "banana"

model[["A"]] <- A
model[["S"]] <- S
model[["F"]] <- F

pList <- omxGenerateParameterList(model)
mList <- omxGenerateSimpleMatrixList(model)
startVals <- omxGenerateValueList(model)
optype <- "cov"
mList[[4]] <- covMatrix
NPSOLOutput <- .Call("callNPSOL", optype, startVals, list(), mList, pList, NA, NA, new.env())
print(NPSOLOutput)
