library(OpenMx)

model <- mxModel("assign",
		 mxMatrix(name="A", nrow=2, ncol=2, free=TRUE, values=1:4),
		 mxMatrix(name="B", nrow=2, ncol=2, free=TRUE, values=0),
		 mxComputeAssign('A', 'B'))
fit <- mxRun(model)
omxCheckEquals(fit@matrices$B@values, fit@matrices$A@values)

model <- mxModel("assign",
                 mxMatrix(name="A", nrow=2, ncol=2, free=c(TRUE,TRUE,FALSE,FALSE), values=1:4),
                 mxMatrix(name="B", nrow=2, ncol=2, free=TRUE, values=0),
                 mxComputeAssign('A', 'B'))
fit <- mxRun(model)
omxCheckEquals(fit@matrices$B@values[,1], fit@matrices$A@values[,1])
omxCheckEquals(fit@matrices$B@values[,2], c(0,0))

model <- mxModel("assign",
                 mxMatrix(name="A", nrow=2, ncol=2, free=TRUE, values=1:4),
                 mxMatrix(name="B", nrow=2, ncol=2, free=FALSE, values=0),
                 mxComputeAssign('A', 'B'))
omxCheckError(mxRun(model), "omxComputeAssign: assign.B[1,1] is fixed, cannot copy assign.A")

model <- mxModel("assign",
                 mxMatrix(name="A", nrow=2, ncol=2, free=TRUE, values=1:4),
                 mxMatrix(name="B", nrow=2, ncol=2, free=FALSE, values=0),
                 mxComputeAssign('A', 'A'))
omxCheckError(mxRun(model), "omxComputeAssign: cannot copy assign.A to itself")

model <- mxModel("assign",
                 mxMatrix(name="A", nrow=2, ncol=2, free=TRUE, values=1:4),
                 mxMatrix(name="B", nrow=2, ncol=2, free=TRUE, values=0),
                 mxMatrix(name="C", nrow=2, ncol=2, free=TRUE, values=0),
                 mxComputeAssign(c('A','A'), c('B','C')))
omxCheckError(mxRun(model), "omxComputeAssign: cannot copy assign.A to more than 1 place")

omxCheckError(mxModel("assign",
                      mxMatrix(name="A", nrow=2, ncol=2, free=TRUE, values=1:4),
                      mxMatrix(name="B", nrow=2, ncol=2, free=TRUE, values=0),
                      mxMatrix(name="C", nrow=2, ncol=2, free=TRUE, values=0),
                      mxComputeAssign(c('A'), c('B','C'))),
              "Arguments 'from' and 'to' must be the same length")
