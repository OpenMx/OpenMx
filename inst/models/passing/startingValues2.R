#
#   Copyright 2007-2018 by the individuals mentioned in the source code history
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.


library(OpenMx)
library(testthat)

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

#-----------------

test_that("mismatched starting values", {
  m2 <- mxModel(mxMatrix("Full", 2, 1, TRUE, values=c(1,2),
                         labels="same", name="A"))
  expect_error(mxRun(m2), "assigned multiple starting values")
  mxRun(omxAssignFirstParameters(m2))
})

test_that("mismatched lbound", {
  m2 <- mxModel(mxMatrix("Full", 2, 1, TRUE, values=1,
                         labels="same", lbound=c(.1,.2), name="A"))
  expect_error(mxRun(m2), "assigned multiple lower bounds")
  mxRun(omxAssignFirstParameters(m2))
})

test_that("mismatched ubound", {
  m2 <- mxModel(mxMatrix("Full", 2, 1, TRUE, values=1,
                         labels="same", ubound=c(.1,.2), name="A"))
  expect_error(mxRun(m2), "assigned multiple upper bounds")
  mxRun(omxAssignFirstParameters(m2))
})

test_that("equal bounds", {
  m2 <- mxModel(mxMatrix("Full", 1, 1, TRUE, values=1.5,
                         lbound=2, ubound=1, name="A"))
  expect_error(mxRun(m2), "greater than or equal to an upper bound")
})

test_that("mxPath bounds", {
  lb <- seq(.75,3,length.out=10)
  lb[6] <- NA
  expect_error(mxPath(letters[1:5], connect = "unique.bivariate",
                      lbound = lb, ubound = 2),
               "2.25 is greater than or equal to upper bound")
})
