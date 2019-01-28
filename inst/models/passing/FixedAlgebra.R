library(OpenMx)

t1 <- mxModel(
  "test",
  mxMatrix("Full", 2,2,values=1:4, free=TRUE, name="z"),
  mxAlgebra(z, name="zCopy", fixed=TRUE, initial=matrix(1., 2,2))
  )

t1 <- mxRun(t1)
omxCheckEquals(t1$zCopy$result, matrix(1., 2,2))

t2 <- mxModel(t1, mxComputeOnce("zCopy", 'fit'))

t2 <- mxRun(t2)
omxCheckEquals(t2$zCopy$result, t2$z$values)

# ---

e1 <- mxModel(
  "test",
  mxMatrix("Full", 2,2,values=1:4, free=TRUE, name="z"),
  mxAlgebra(z, name="zCopy", fixed=TRUE))

omxCheckError(mxRun(e1), "Missing value found in initial result of algebra 'test.zCopy'")

e2 <- mxModel(
  "test",
  mxMatrix("Full", 2,2,values=1:4, free=TRUE, name="z"),
  mxAlgebra(z, name="zCopy", initial=matrix(pi,1,1)))

omxCheckError(mxRun(e2), "In algebra 'test.zCopy' initial result will be immediately discarded")
