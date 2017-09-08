library(OpenMx)

df <- data.frame(dv=c(1L:10L,NA))

m1 <- mxModel(
  "integer",
  mxData(df, 'raw'),
  mxMatrix(nrow=1,ncol=1,labels="data.dv",name="obj"),
  mxFitFunctionAlgebra("obj"))

omxCheckError(mxRun(m1), "integer.data: NA in definition variable 'dv' row 11")

df <- data.frame(dv=c(2.5, 3.5, NA))

m1 <- mxModel(
  "double",
  mxData(df, 'raw'),
  mxMatrix(nrow=1,ncol=1,labels="data.dv",name="obj"),
  mxFitFunctionAlgebra("obj"))

omxCheckError(mxRun(m1), "double.data: NA in definition variable 'dv' row 3")
