library(OpenMx)

testFit <- function(objective, startVals=c(), bounds=c(), matList=list(), varList=list(), data=c(), state=c()) {
	return(.Call("callNPSOL", objective, startVals, bounds, matList, varList, data, state));
}