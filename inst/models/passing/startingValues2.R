library(OpenMx)

# Try with 2 appearances of 'x', one starting at 1.0 and the other NA

m1 <- mxModel(mxMatrix(nrow=1,ncol=1,free=TRUE,labels="x",name="mat1"),
              mxMatrix(nrow=1,ncol=1,free=TRUE,values=1,labels="x",name="mat2"),
              mxAlgebra(abs(mat1), "abs"),
              mxFitFunctionAlgebra("abs"))

m1$mat1$values[1,1] <- NA

afp <- omxAssignFirstParameters(m1)
omxCheckEquals(afp$mat1$values[1,1], 1)
omxCheckEquals(afp$mat2$values[1,1], 1)

params <- omxGetParameters(m1)
omxCheckEquals(params, 1)

m1 <- mxRun(m1)
omxCheckCloseEnough(m1$output$estimate, 0, 1e-6)
omxCheckEquals(m1$mat1$values, m1$mat2$values)

#----------------- single appearance of 'x', no starting value

m1 <- mxModel(mxMatrix(nrow=1,ncol=1,free=TRUE,labels="x",name="mat1"),
              mxAlgebra(abs(mat1), "abs"),
              mxFitFunctionAlgebra("abs"))

m1$mat1$values[1,1] <- NA
omxCheckError(mxRun(m1), "Parameter 'x' has no starting value")

#----------------- two appearances of 'x', no starting values

m1 <- mxModel(mxMatrix(nrow=1,ncol=1,free=TRUE,labels="x",name="mat1"),
              mxMatrix(nrow=1,ncol=1,free=TRUE,labels="x",name="mat2"),
              mxAlgebra(abs(mat1), "abs"),
              mxFitFunctionAlgebra("abs"))

m1$mat1$values[1,1] <- NA
m1$mat2$values[1,1] <- NA

omxCheckError(mxRun(m1), "Parameter 'x' has no starting value")
