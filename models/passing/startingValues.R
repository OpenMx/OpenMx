library(OpenMx)

covMatrix <- matrix( c(0.77642931, 0.39590663, 0.39590663, 0.49115615), 
	nrow = 2, ncol = 2, byrow = TRUE, dimnames = list(c('a','b'), c('a','b')))

model <- mxModel("missingtest",
	   mxMatrix("Full", values = c(0,0.2,0,0), name="A", nrow=2, ncol=2),
	   mxMatrix("Symm", values = c(0.8,0,0,0.8), name="S", nrow=2, ncol=2, free=TRUE),
	   mxMatrix("Iden", name="F", nrow=2, ncol=2, dimnames = list(c('a','b'), c('a','b'))),
	   mxData(covMatrix, 'cov', numObs = 100),
	   mxExpectationRAM("A", "S", "F"),
	   mxFitFunctionML())

model[["A"]]@free[2,1] <- TRUE
model[["A"]]@values[2,1] <- NA   # oops

msg <<- ''
tryCatch(mxRun(model, silent=TRUE), error = function(e) msg <<- e$message)
omxCheckEquals(msg, "Starting value in missingtest.A[2,1] is missing")
