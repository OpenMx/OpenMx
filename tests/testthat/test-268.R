library(testthat)
library(OpenMx)
data(demoOneFactor)
context("268")

body <- function() {
  E <- mxAlgebra(U%x%1,name="E")
  
  factorModel <- mxModel(
    "One Factor",
    mxMatrix("Full", 5, 1, values=0.8,
             free=TRUE, name="A"),
    mxMatrix("Symm", 1, 1, values=1,
             free=FALSE, name="L"),
    mxMatrix("Diag", 5, 5, values=1,
             free=TRUE, name="U"),
    mxAlgebra(A %*% L %*% t(A) + E, name="R"),
    mxExpectationNormal(covariance = "R",
                        dimnames = names(demoOneFactor)),
    mxFitFunctionML(),
    mxData(cov(demoOneFactor), type="cov", numObs=500)
  )
  
  expect_error(mxCheckIdentification(factorModel,details=T),
               "'E' not found")
}

test_that("268", body)
