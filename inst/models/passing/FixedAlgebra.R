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
