library(testthat)
context("discrete")
library(OpenMx)

# to test:
# Normal, RAM, LISREL
# with and without regular thresholds

skip("not yet")

factorModel <- mxModel(
  "One Factor",
  mxMatrix("Full", 5, 1, values=0.8,
           free=TRUE, name="A"),
  mxMatrix("Symm", 1, 1, values=1,
           free=FALSE, name="L"),
  mxMatrix("Diag", 5, 5, values=1,
           free=TRUE, name="U"),
  mxAlgebra(A %*% L %*% t(A) + U, name="R"),
  mxMatrix(nrow=5, ncol=3,
           values=c(4, -2, 1, .7, NA,
                    5, -2, 1, .6, NA,
                    6, -2, 2, 4, .5),
           name="D"),
  mxExpectationNormal(covariance = "R",
                      dimnames = paste0('x',1:5),
                      threshnames = paste0('x',1:3),
                      discrete = "D"),
  mxFitFunctionML())

mxGetExpected(factorModel, "thresholds")

factorModel <- mxGenerateData(factorModel, 200, returnModel = TRUE)

# convert to raw counts
factorModel$data$observed$x1 <-
  unclass(factorModel$data$observed$x1) - 1L

fit <- mxRun(factorModel)

stop("here")

# ---------------

library(OpenMx)

manifests <- paste0('x',1:5)
latents <- c("G")
factorModel <- mxModel(
  "One Factor",
  type="RAM",
  manifestVars = manifests,
  latentVars = latents,
  mxPath(from=latents, to=manifests,values=0.8),
  mxPath(from=manifests, arrows=2,values=1),
  mxPath(from=latents, arrows=2,
         free=FALSE, values=1.0),
  mxPoisson(paste0("x",1:2), c(4,5), c(.6,.7)),
  mxNegativeBinomial("x3", 6, 4, .5))

round(mxGetExpected(factorModel, "thresholds"),3)

mxGenerateData(factorModel, 20)
