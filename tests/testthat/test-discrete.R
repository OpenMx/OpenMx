library(testthat)
context("discrete")
library(OpenMx)

# to test:
# Normal, RAM, LISREL
# ML / WLS
# with and without regular thresholds
# different column orders for discrete vs discreteSpec
# check equivalence of parameterizations (total var = 1 vs loading = 1)
# ensure front and back-end have exactly the same expectations
# focus on mxMarginalNegativeBinomial recovery
# ordered factor vs raw count

test_that("Normal", {
  library(OpenMx)
  library(testthat)
  
  RNGversion("4.0")
  set.seed(1)
  
  factorModel <- mxModel(
    "One Factor",
    mxMatrix(nrow=1, ncol=5, free=FALSE, values=0, name="M"),
    mxMatrix("Full", 5, 1, values=0.8,
             free=TRUE, name="A"),
    mxMatrix("Symm", 1, 1, values=1,
             free=FALSE, name="L"),
    mxMatrix("Diag", 5, 5, values=1,
             free=FALSE, name="U"),
    mxAlgebra(A %*% L %*% t(A) + U, name="R"),
    mxMatrix(nrow=3, ncol=3,
             values=c(.01, 1.1, NA,
                      .01, 2, NA,
                      .01, 4, .5),
             dimnames=list(c(), paste0('x',1:3)),
             name="D"),
    mxExpectationNormal(covariance = "R",
                        dimnames = paste0('x',1:5),
                        discrete = "D",
                        discreteSpec = matrix(c(4, 1, 5, 1, 6, 2), ncol=3,
                                              dimnames=list(c(), paste0('x',1:3))),
                        means = "M"),
    mxFitFunctionML())
  
  thr <- mxGetExpected(factorModel, "thresholds")
  expect_equal(thr[1:4,'x1'], c(-0.53, 0.68, 1.65, 2.5), .01)
  expect_equal(thr[1:5,'x2'], c(-1.36, -0.28, 0.6, 1.38, 2.08), .01)
  expect_equal(thr[1:6,'x3'], c(-1.87, -1.1, -0.49, 0.02, 0.46, 0.86), .01)
  
  factorModel <- mxGenerateData(factorModel, 400, returnModel = TRUE)
  
  # convert to raw counts
  # factorModel$data$observed$x1 <-
  #   unclass(factorModel$data$observed$x1) - 1L
  
  factorModel$D$free <- !is.na(factorModel$D$values)
  
  fit <- mxRun(factorModel)
  expect_equal(fit$output$fit, 6256.76, .01)
  dv <- fit$D$values
  dv <- dv[!is.na(dv)]
  expect_equal(dv, c(0, 1.047, 0.022, 2.072, 0.013, 2.91, 0.411), .01)
})

# ---------------

test_that("RAM", {
  library(OpenMx)
  library(testthat)
  
  RNGversion("4.0")
  set.seed(1)
  
  manifests <- paste0('x',1:5)
  latents <- c("G")
  factorModel <- mxModel(
    "One Factor",
    type="RAM",
    manifestVars = manifests,
    latentVars = latents,
    mxPath(from=latents, to=manifests,values=0.8),
    mxPath(from=manifests, arrows=2,values=1, free=FALSE),
    mxPath(from=latents, arrows=2,
           free=FALSE, values=1.0),
    mxPath(from = 'one', to = manifests, free=FALSE),
    mxMarginalPoisson(paste0("x",1:2), c(4,5), c(1.1, 2)),
    mxMarginalNegativeBinomial("x3", 6, 4, .5))
  
  trueDv <- factorModel$Discrete$values
  trueDv <- trueDv[!is.na(trueDv)]
  
  thr <- mxGetExpected(factorModel, "thresholds")
  expect_equal(thr[1:4,'x1'], c(-0.53, 0.68, 1.65, 2.5), .01)
  expect_equal(thr[1:5,'x2'], c(-1.36, -0.28, 0.6, 1.38, 2.08), .01)
  expect_equal(thr[1:6,'x3'], c(-1.87, -1.1, -0.49, 0.02, 0.46, 0.86), .01)
  
  factorModel <- mxGenerateData(factorModel, 400, returnModel = TRUE)
  
  # autodetect maximum count from data
  factorModel$expectation$discreteSpec[1,] <- NA
  
  if (0) {
    m1 <- mxRun(mxModel(factorModel, mxComputeSequence(list(
      mxComputeOnce('expectation'),
      mxComputeReportExpectation()))))
    
    m1$expectation$output
  }
  
  factorModel$Discrete$free <- !is.na(factorModel$Discrete$values)
  
  fit <- mxRun(factorModel)
  #  summary(fit)
  expect_equal(fit$output$fit, 6256.76, .01)
  dv <- fit$Discrete$values
  dv <- dv[!is.na(dv)]
  expect_equal(dv, c(0, 1.047, 0.022, 2.072, 0.013, 2.91, 0.411), .01)
})
